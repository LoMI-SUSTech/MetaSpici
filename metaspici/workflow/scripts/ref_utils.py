#!/usr/bin/env python3
"""Validated input conversion and pair reconstruction helpers."""
import argparse
from collections import Counter
import gzip
import itertools
import json
from pathlib import Path
import re
import sqlite3
import sys

def fasta(path):
    name, seq = None, []
    with open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq).upper()
                name, seq = line[1:].split()[0], []
            elif line.strip():
                if name is None:
                    raise ValueError(f'Sequence before FASTA header: {path}')
                seq.append(''.join(line.split()))
    if name is not None:
        yield name, ''.join(seq).upper()

def validate_fasta(path):
    seen = set()
    for name, seq in fasta(path):
        if name in seen or not seq or not re.fullmatch('[ACGTRYSWKMBDHVN]+', seq):
            raise ValueError(f'Duplicate ID, empty or invalid DNA: {path}: {name}')
        seen.add(name)
    if not seen:
        raise ValueError(f'Empty FASTA: {path}')

def validate_inputs(mags, mge):
    genome_names = set()
    for path in mags:
        validate_fasta(path)
        genome_names.add(Path(path).stem)
    validate_fasta(mge)
    overlap = genome_names & {name for name, seq in fasta(mge)}
    if overlap:
        raise ValueError('MAG and MGE IDs collide: ' + ', '.join(sorted(overlap)))

def merge_genomes(directory, output, offsets_path=None):
    paths = sorted(Path(directory).glob('*.fa'))
    if not paths:
        raise ValueError('No genomes passed QC/dereplication')
    with open(output, 'w') as handle, open(offsets_path or output + '.contig_offsets.tsv', 'w') as offsets:
        offsets.write('genome\tcontig\tstart0\tend0\n')
        for path in paths:
            # Preserve the original genome-level concatenation model, with documented offsets.
            # No padding is introduced here (which would change abundance denominators).
            records = list(fasta(path))
            start = 0
            for name, seq in records:
                offsets.write(f'{path.stem}\t{name}\t{start}\t{start+len(seq)}\n')
                start += len(seq)
            handle.write('>' + path.stem + '\n' + ''.join(s for _, s in records) + '\n')

def fastq(path):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt') as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq, plus, qual = (handle.readline().rstrip('\r\n') for _ in range(3))
            if not header.startswith('@') or not plus.startswith('+') or len(seq) != len(qual) or not seq:
                raise ValueError(f'Malformed FASTQ: {path}')
            name = re.sub(r'/[12]$', '', header.split()[0][1:])
            yield name, seq, qual

def normalize(r1, r2, out1, out2):
    count = 0
    with open(out1, 'w') as x, open(out2, 'w') as y:
        for a, b in itertools.zip_longest(fastq(r1), fastq(r2)):
            if a is None or b is None or a[0] != b[0]:
                raise ValueError(f'Paired FASTQ IDs/counts differ at pair {count+1}')
            for out, rec in ((x, a), (y, b)):
                out.write(f'@{rec[0]}\n{rec[1]}\n+\n{rec[2]}\n')
            count += 1
    if not count:
        raise ValueError('No reads survived fastp')

def mask_bed(mode, source, dest):
    with open(source) as handle, open(dest, 'w') as out:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip().split('\t')
            try:
                if mode == 'gff':
                    name, start, end = parts[0], int(parts[3]), int(parts[4])
                elif mode == 'metaxa':
                    name, start, end = parts[0], int(parts[8]), int(parts[9])
                elif mode == 'dust':
                    # NCBI dustmasker -outfmt interval: 0-based inclusive.
                    if line.startswith('>'):
                        current = line[1:].split()[0]
                        continue
                    start, end = [int(x.strip()) for x in line.split('-')]
                    out.write(f'{current}\t{start}\t{end+1}\n')
                    continue
                else:
                    raise ValueError(mode)
            except (IndexError, ValueError) as e:
                raise ValueError(f'Cannot parse {mode} row: {line.rstrip()}') from e
            start, end = sorted((start, end))
            if start < 1:
                raise ValueError(f'Invalid 1-based {mode} coordinate: {line.rstrip()}')
            out.write(f'{name}\t{start-1}\t{end}\n')

def representatives(fna, clusters, output):
    selected = set()
    with open(clusters) as handle:
        for line in handle:
            if line.strip():
                selected.add(line.split('\t')[0])
    found = set()
    with open(output, 'w') as out:
        for name, seq in fasta(fna):
            if name in selected:
                found.add(name)
                out.write(f'>{name}\n{seq}\n')
    if not found or found != selected:
        raise ValueError('Invalid/empty viral clustering representatives')

def contacts(source, output, audit):
    """Input BED sorted by normalized full QNAME; retain pairs, then coordinate dedup.

    Uses SQLite for duplicate coordinates, avoiding an unbounded in-memory set.
    Two distinct primary mates must survive; counts are unordered genome pairs.
    """
    db_path = output + '.dedup.sqlite'
    Path(db_path).unlink(missing_ok=True)
    db = sqlite3.connect(db_path)
    db.execute('CREATE TABLE seen (key TEXT PRIMARY KEY) WITHOUT ROWID')
    counts, metrics = Counter(), Counter()
    def name(line):
        return re.sub(r'/[12]$', '', line.split('\t')[3])
    with open(source) as handle:
        previous = None
        for read, lines in itertools.groupby(handle, name):
            if previous is not None and read < previous:
                raise ValueError('BED must be sorted by normalized QNAME')
            previous = read
            rows = [x.rstrip().split('\t') for x in lines]
            metrics['read_ids_seen'] += 1
            if len(rows) != 2:
                metrics['not_exactly_two_alignments'] += 1
                continue
            # bedtools bamtobed appends /1 and /2 to paired primary alignments.
            mates = {r[3][-2:] for r in rows}
            if mates != {'/1', '/2'}:
                raise ValueError('Expected one /1 and one /2 BED alignment for ' + read)
            coords = sorted((r[0], int(r[1])) for r in rows)
            key = json.dumps(coords, separators=(',', ':'))
            before = db.total_changes
            db.execute('INSERT OR IGNORE INTO seen VALUES (?)', (key,))
            if db.total_changes == before:
                metrics['coordinate_duplicate_pairs'] += 1
                continue
            counts[tuple(sorted(r[0] for r in rows))] += 1
            metrics['retained_pairs'] += 1
    db.commit()
    db.close()
    Path(db_path).unlink()
    with open(output, 'w') as out:
        for (a, b), n in sorted(counts.items()):
            out.write(f'{a}\t{b}\t{n}\n')
    Path(audit).write_text(json.dumps(dict(metrics), indent=2))

if __name__ == '__main__':
    command, *args = sys.argv[1:]
    {'normalize': normalize, 'mask_bed': mask_bed,
     'merge_genomes': merge_genomes, 'representatives': representatives,
     'contacts': contacts}[command](*args)
