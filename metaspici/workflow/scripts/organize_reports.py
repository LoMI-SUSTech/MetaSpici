"""Collect small tool reports and inventory the completed run without moving intermediates."""
import argparse
import csv
from pathlib import Path
import shutil


def organize(root, sample):
    root = Path(root)
    reports = root / 'qc/reference/tool_reports'
    hicup = root / 'qc/reads/hicup'
    reports.mkdir(parents=True, exist_ok=True)
    hicup.mkdir(parents=True, exist_ok=True)
    (root / 'run_info').mkdir(parents=True, exist_ok=True)
    # Keep native tool outputs intact; only copy small, useful reports for convenient review.
    for name in ('Bdb.csv', 'Chdb.csv', 'Cdb.csv', 'Wdb.csv'):
        src = root / 'intermediates/reference_qc/drep/data_tables' / name
        if not src.is_file():
            continue
        dest = reports / 'drep' / src.name
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
    for name in ('quality_summary.tsv', 'completeness.tsv', 'contamination.tsv', 'complete_genomes.tsv'):
        src = root / 'intermediates/reference_qc/checkv' / name
        if src.is_file():
            dest = reports / 'checkv' / name
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dest)
    for src in sorted((root / 'intermediates/reads/truncated').glob('*')):
        if src.is_file() and src.name.lower().endswith(('.txt', '.tsv', '.html', '.json')):
            shutil.copy2(src, hicup / src.name)
    guide = f'''# MetaSpici outputs: {sample}

One output directory represents one sample/run. Parameters are recorded in `run_info/config.yaml`.

| Location | Contents |
|---|---|
| `results/{sample}.associations.unfiltered_fe.tsv` | All non-self associations after depth filtering; FE values retained; no FE/contact-count cutoff |
| `results/{sample}.mge_host.filtered.tsv` | MGE-host associations passing FE and minimum contacts; source=MGE, target=host |
| `qc/network/{sample}.genome_depth.tsv` | Per-genome read-end counts and depth before depth filtering |
| `qc/network/pair_filtering.json` | Pair reconstruction and coordinate deduplication audit |
| `qc/reads/fastp.html`, `fastp.json` | Read QC reports |
| `qc/reads/hicup/` | Native HiCUP summary reports (names may include timestamps) |
| `qc/reference/mge_clusters.tsv` | MGE clustering assignments, present in the --votu branch only |
| `qc/reference/tool_reports/drep/` | Copies of selected native dRep data tables (Bdb/Chdb/Cdb/Wdb) |
| `qc/reference/tool_reports/checkv/` | Copies of CheckV reports, present in the --votu branch only |
| `reference/host_genomes.fasta` | Retained hosts, each genome concatenated as before |
| `reference/host_contig_offsets.tsv` | Original contig coordinates within concatenated hosts |
| `reference/mge_sequences.fasta` | Retained MGE sequences |
| `reference/combined.fasta` | Combined host and MGE reference |
| `reference/combined.mask1.fasta` | Reference after Mask1 |
| `reference/masking/mask1_regions.bed`, `mask2_regions.bed` | Masking/filtering intervals |
| `intermediates/reads/` | fastp, normalized, fastuniq and truncated FASTQ files |
| `intermediates/reference_qc/` | Complete native dRep/CheckV directories and MGE BLAST/ANI work files |
| `intermediates/masking/` | Metaxa, rRNA/tRNA annotations, low-complexity files and Mask2 BLAST data |
| `intermediates/alignment/` | BWA index, BAM and pre-Mask2 BED files |
| `intermediates/network/alignments.mask2_filtered.qname_sorted.bed` | Mask2-filtered, read-name-sorted alignments |
| `intermediates/network/contact_counts.tsv` | Undirected contact counts including self contacts; columns source, target, count (no header) |
| `logs/` | Per-step stdout/stderr logs; stdout-as-data steps log stderr separately |
| `run_info/config.yaml` | Effective input paths and analysis parameters |
| `run_info/invocation.json` | Invocation arguments and Snakemake command |
| `run_info/output_manifest.tsv` | Snapshot inventory of actual files, excluding itself and workflow engine bookkeeping |
| `assembly_workflow/` | Numbered assembly-only stages: `01_assembly/`, `02_binning/`, `03_mge_detection/`, `04_hic_network/`, `05_mge_host_linking/`; plus assembly-specific `qc/`, `logs/`, and `run_info/` |

Default thresholds: fixed Depth_Total >=0.15, FE >50, final contact count >=4; actual values are in config.yaml.
Self contacts contribute to FE calculations but are absent from both association tables.
Tool-native filenames are preserved inside their own directories; large intermediates are not duplicated or deleted.
The manifest is a final-report snapshot, not a checksum-based integrity certificate. Logs may still receive final lines afterward.
'''
    (root / 'README_outputs.md').write_text(guide)
    categories = {'results':'association_results', 'qc':'quality_reports', 'reference':'reference_sequences_and_intervals',
                  'intermediates':'intermediate_data', 'logs':'step_logs', 'run_info':'run_provenance',
                  'assembly_workflow':'assembly_workflow_artifacts'}
    manifest = root / 'run_info/output_manifest.tsv'
    with manifest.open('w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['relative_path', 'category', 'size_bytes'])
        for path in sorted(root.rglob('*')):
            rel = path.relative_to(root)
            if not path.is_file() or path == manifest or any(part.startswith('.') for part in rel.parts):
                continue
            writer.writerow([rel.as_posix(), categories.get(rel.parts[0], 'documentation'), path.stat().st_size])


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True)
    p.add_argument('--root', type=Path, default=Path('.'))
    a = p.parse_args()
    organize(a.root, a.sample)
