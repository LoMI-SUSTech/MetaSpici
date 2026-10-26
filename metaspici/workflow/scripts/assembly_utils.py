#!/usr/bin/env python3
"""Small deterministic helpers for the short-read + Hi-C assembly workflow."""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

DNA = re.compile(r"^[ACGTRYSWKMBDHVN]+$", re.IGNORECASE)


def read_fasta(path: Path):
    name = None
    seq = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            else:
                if name is None:
                    raise ValueError(f"Sequence occurs before a FASTA header in {path}")
                seq.append("".join(line.split()))
    if name is not None:
        yield name, "".join(seq).upper()


def validate_fasta(path: Path) -> int:
    seen = set()
    count = 0
    for name, seq in read_fasta(path):
        if name in seen:
            raise ValueError(f"Duplicate FASTA ID {name!r} in {path}")
        if not seq or not DNA.fullmatch(seq):
            raise ValueError(f"Empty or invalid DNA sequence {name!r} in {path}")
        seen.add(name)
        count += 1
    if count == 0:
        raise ValueError(f"FASTA is empty: {path}")
    return count


def fasta_files(directory: Path, recursive: bool = False):
    if not directory.is_dir():
        raise ValueError(f"Bin directory does not exist: {directory}")
    candidates = directory.rglob("*") if recursive else directory.iterdir()
    files = sorted(p for p in candidates if p.is_file() and p.suffix.lower() in {".fa", ".fna", ".fasta"})
    return files


def combine_bins(hic_dir: Path, short_dir: Path, output_dir: Path) -> int:
    """Prefix bin filenames by source to avoid collisions; preserve FASTA records."""
    hic = fasta_files(hic_dir, recursive=True)
    short = fasta_files(short_dir, recursive=True)
    if not hic and not short:
        raise ValueError("Neither MetaTOR nor BASALT produced any FASTA bins")
    output_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for source, paths in (("hic", hic), ("short", short)):
        for path in paths:
            stem = path.name
            for suffix in (".fasta", ".fna", ".fa"):
                if stem.lower().endswith(suffix):
                    stem = stem[:-len(suffix)]
                    break
            dest = output_dir / f"{source}_{stem}.fa"
            if dest.exists():
                raise ValueError(f"Refusing to overwrite duplicate combined bin name: {dest.name}")
            shutil.copy2(path, dest)
            validate_fasta(dest)
            count += 1
    return count


def resolve_basalt_binset(workdir: Path) -> Path:
    """Resolve a known final-binset directory name across the legacy BASALT output contract."""
    workdir = workdir.resolve()
    preferred = (
        "BestBinset_outlier_refined_filtrated_retrieved_retrieved",
        "BestBinset_outlier_refined_filtrated_retrieved",
        "BestBinset_outlier_refined",
        "Final_bestbinset",
        "Final_BestBinset",
        "BestBinset",
        "Final_binset",
    )
    for name in preferred:
        candidate = workdir / name
        if candidate.is_dir() and fasta_files(candidate):
            return candidate
    raise ValueError(
        "BASALT finished without a recognized final binset directory. Expected one of: "
        + ", ".join(preferred)
        + f". Directories found: {[p.name for p in workdir.iterdir() if p.is_dir()]}"
    )


def validate_generated_inputs(mag_dir: Path, mge_fasta: Path):
    mags = fasta_files(mag_dir)
    if not mags:
        raise ValueError(f"No MAG FASTA files were produced in {mag_dir}")
    names = {p.stem for p in mags}
    for p in mags:
        validate_fasta(p)
    mge_count = validate_fasta(mge_fasta)
    mge_ids = {name for name, _ in read_fasta(mge_fasta)}
    overlap = names & mge_ids
    if overlap:
        raise ValueError("MAG filenames collide with MGE FASTA IDs: " + ", ".join(sorted(overlap)))
    return len(mags), mge_count



def combine_bins_many(output_dir: Path, sources: dict[str, Path]) -> int:
    """Copy multiple binsets using a stable prefix per source (including 1 Mb contigs)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    count = 0
    for label, directory in sources.items():
        files = fasta_files(directory, recursive=True)
        for index, path in enumerate(files, start=1):
            dest = output_dir / f'{label}_{index:06d}.fa'
            shutil.copy2(path, dest)
            validate_fasta(dest)
            count += 1
    if not count:
        raise ValueError('No MAG candidates found across all supplied bin sources')
    return count


def split_long_contigs(fasta: Path, outdir: Path, min_length: int = 1000000) -> int:
    """Make each >=1 Mb contig a separate MAG candidate; zero is permitted."""
    outdir.mkdir(parents=True, exist_ok=True)
    count = 0
    for name, seq in read_fasta(fasta):
        if len(seq) < min_length:
            continue
        count += 1
        (outdir / f'long_{count:06d}.fa').write_text(f'>{name}\n{seq}\n')
    return count


def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    combine = sub.add_parser("combine-bins")
    combine.add_argument("hic_dir", type=Path)
    combine.add_argument("short_dir", type=Path)
    combine.add_argument("output_dir", type=Path)
    resolve = sub.add_parser("resolve-basalt-binset")
    resolve.add_argument("workdir", type=Path)
    combine_many = sub.add_parser("combine-bins-many")
    combine_many.add_argument("output_dir", type=Path)
    combine_many.add_argument("sources", nargs='+', help='label=directory (can be empty for absent branch)')
    split_large = sub.add_parser("split-long-contigs")
    split_large.add_argument("fasta", type=Path)
    split_large.add_argument("output_dir", type=Path)
    split_large.add_argument("--min-length", type=int, default=1000000)
    validate = sub.add_parser("validate-inputs")
    validate.add_argument("mag_dir", type=Path)
    validate.add_argument("mge_fasta", type=Path)
    args = parser.parse_args()
    if args.command == "combine-bins":
        print(f"Copied {combine_bins(args.hic_dir, args.short_dir, args.output_dir)} candidate bins")
    elif args.command == "combine-bins-many":
        sources = {}
        for source in args.sources:
            if '=' not in source:
                raise ValueError('Expected source in label=directory format')
            label, directory = source.split('=', 1)
            if directory:
                sources[label] = Path(directory)
        print(f'Copied {combine_bins_many(args.output_dir, sources)} MAG candidates')
    elif args.command == "split-long-contigs":
        print(f'Extracted {split_long_contigs(args.fasta, args.output_dir, args.min_length)} contigs')
    elif args.command == "resolve-basalt-binset":
        print(resolve_basalt_binset(args.workdir))
    else:
        n_mags, n_mges = validate_generated_inputs(args.mag_dir, args.mge_fasta)
        print(f"Validated {n_mags} MAG FASTA files and {n_mges} MGE sequences")


if __name__ == "__main__":
    main()
