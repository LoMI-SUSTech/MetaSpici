#!/usr/bin/env python3
"""Run the selected assembly/MGE stage, then the shared reference association workflow."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import yaml

ROOT = Path(__file__).resolve().parent
SNAKEFILE = ROOT / 'Snakefile_assembly'
SNAKEFILE_EXTENDED = ROOT / 'Snakefile_assembly_extended'
LONG_TYPES = ('nanopore', 'pacbio', 'nanoporeR10', 'pacbioHiFi')
REF_SCRIPT = ROOT / 'metaspici_ref.py'


def conda_environment_variables(env_root, inherited=None):
    env = dict(os.environ if inherited is None else inherited)
    root = str(Path(env_root).expanduser().resolve())
    paths = [root]
    for item in env.get('CONDA_ENVS_PATH', '').split(os.pathsep):
        item = item.strip()
        if not item:
            continue
        resolved = str(Path(item).expanduser().resolve())
        if resolved not in paths:
            paths.append(resolved)
    env['CONDA_ENVS_PATH'] = os.pathsep.join(paths)
    return env


def check_inputs(args):
    required_reads = ['hic_r1', 'hic_r2']
    if args.assembly_variant in ('short_hic', 'hybrid_hic'):
        required_reads += ['short_r1', 'short_r2']
    if args.assembly_variant in ('long_hic', 'hybrid_hic'):
        required_reads += ['long_reads']
        if args.long_type not in LONG_TYPES:
            raise ValueError('--long-type must be one of: ' + ', '.join(LONG_TYPES))
    if args.assembly_variant == 'hybrid_hic':
        opera = args.opera_ms
        if not opera:
            raise ValueError('Hybrid mode requires a complete OPERA-MS installation. '
                             'Provide --opera-ms /absolute/path/to/OPERA-MS.pl; '
                             'the original standalone Perl file cannot work without its upstream bin/ and tools.')
        opera = Path(opera).expanduser().resolve()
        if not opera.is_file() or not (opera.parent / 'bin').is_dir():
            raise ValueError(f'OPERA-MS installation missing launcher or sibling bin/ directory: {opera}')
        args.opera_ms = str(opera)
    for label in required_reads:
        path = Path(getattr(args, label)).expanduser().resolve()
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f'{label.replace("_", "-")} is missing or empty: {path}')
        setattr(args, label, str(path))
    for label in ('checkm_db', 'checkv_db', 'genomad_db'):
        path = Path(getattr(args, label)).expanduser().resolve()
        if not path.is_dir():
            raise ValueError(f'{label.replace("_", "-")} directory not found: {path}')
        setattr(args, label, str(path))
    if not (Path(args.checkm_db) / 'genome_tree/genome_tree.derep.txt').is_file():
        raise ValueError(f'CheckM database is incomplete: {args.checkm_db}')
    if not (Path(args.basalt_prefix) / 'bin/BASALT').is_file():
        raise ValueError(f'BASALT executable not found: {args.basalt_prefix}/bin/BASALT')
    if args.threads < 1:
        raise ValueError('--threads must be a positive integer')
    if not 95 <= args.identity <= 99:
        raise ValueError('--identity must be between 95 and 99 (percent)')
    if not args.sample or not __import__('re').fullmatch(r'[A-Za-z0-9_-]+', args.sample):
        raise ValueError('--sample must contain only letters, digits, underscore or hyphen')


def atomic_yaml(path: Path, cfg: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.yaml.tmp')
    tmp.write_text(yaml.safe_dump(cfg, sort_keys=False))
    tmp.replace(path)


def prepare_config(args, out: Path):
    legacy_stage = out / 'assembly_stage'
    stage = out / 'assembly_workflow'
    cfg_path = stage / 'config.yaml'
    if (legacy_stage / 'config.yaml').exists() and not cfg_path.exists():
        raise ValueError('This output directory contains a legacy assembly_stage layout. Use a new --outdir to avoid mixing output layouts.')
    input_signatures = {}
    for key in ('short_r1', 'short_r2', 'hic_r1', 'hic_r2', 'long_reads'):
        if getattr(args, key, None) is None:
            continue
        item = Path(getattr(args, key)).stat()
        input_signatures[key] = {'size_bytes': item.st_size, 'mtime_ns': item.st_mtime_ns}
    cfg = {
        'root': str(ROOT),
        'sample': args.sample,
        'input_signatures': input_signatures,
        'threads': args.threads,
        'short_read1': args.short_r1,
        'short_read2': args.short_r2,
        'hic_read1': args.hic_r1,
        'hic_read2': args.hic_r2,
        'enzyme': args.enzyme,
        'identity': args.identity,
        'checkm_db': args.checkm_db,
        'checkv_db': args.checkv_db,
        'genomad_db': args.genomad_db,
        'basalt_prefix': str(Path(args.basalt_prefix).expanduser().resolve()),
        'basalt_install_source_ref': 'v1.1.0 (used only when setup installs BASALT)',
        'source_cache_dir': str(Path(args.cache_dir).expanduser().resolve()),
        'conda_executable': str(Path(args.conda_executable).expanduser().resolve()),
        'output_root': str(out),
        'software_version': '1.1',
    }
    # Keep the established Short + Hi-C workflow configuration contract.
    if args.assembly_variant != 'short_hic':
        cfg.update({
            'assembly_variant': args.assembly_variant,
            'long_reads': args.long_reads,
            'long_type': args.long_type,
            'opera_ms': args.opera_ms,
            'software_version': '1.1',
        })
    active_snakefile = SNAKEFILE if args.assembly_variant == 'short_hic' else SNAKEFILE_EXTENDED
    if cfg_path.exists():
        old = yaml.safe_load(cfg_path.read_text()) or {}
        ignored = {'threads', 'dry_run'}
        changed = [k for k, v in cfg.items() if k not in ignored and old.get(k) != v]
        if changed:
            raise ValueError('Assembly output directory already contains a run with different settings (' + ', '.join(changed) + '); choose a new --outdir')
    atomic_yaml(cfg_path, cfg)
    invocation = {
        'argv': sys.argv,
        'snakemake': ['snakemake', '-s', str(active_snakefile), '--directory', str(stage),
                      '--configfile', str(cfg_path), '--cores', str(args.threads),
                      '--use-conda', '--printshellcmds', '--rerun-incomplete'],
        'reference_args': args.reference_args,
    }
    info = stage / 'run_info'
    info.mkdir(parents=True, exist_ok=True)
    (info / 'invocation.json').write_text(json.dumps(invocation, indent=2))
    return stage, cfg_path, cfg


def run_snakemake(stage: Path, cfg_path: Path, args, env):
    snakefile = SNAKEFILE if args.assembly_variant == 'short_hic' else SNAKEFILE_EXTENDED
    command = ['snakemake', '-s', str(snakefile), '--directory', str(stage),
               '--configfile', str(cfg_path), '--cores', str(args.threads),
               '--use-conda', '--printshellcmds', '--rerun-incomplete']
    if args.dry_run:
        command.append('--dry-run')
    print('Assembly-stage command:', ' '.join(command), flush=True)
    subprocess.run(command, check=True, env=env)


def ref_command(args, out: Path, mag_dir: Path, mge_fasta: Path, extra, dry_run=False):
    command = [sys.executable, str(REF_SCRIPT),
               '--env-root', str(Path(args.env_root).resolve()),
               '--mag-dir', str(mag_dir.resolve()),
               '--mge-fasta', str(mge_fasta.resolve()),
               '--hic-read1', args.hic_r1, '--hic-read2', args.hic_r2,
               '--sample', args.sample, '--threads', str(args.threads),
               '--identity', f'{args.identity:g}', '--re1', f'^GATC,{args.enzyme}',
               '--checkm-db', args.checkm_db, '--checkv-db', args.checkv_db,
               '--output', str(out.resolve())]
    if dry_run:
        command.append('--dry-run')
    command.extend(extra)
    return command


def reference_dry_run(args, out: Path, extra, env):
    # Assembly products do not exist yet during a dry run. Use tiny valid placeholder
    # FASTAs in a temporary directory to independently check the shared reference workflow DAG.
    with tempfile.TemporaryDirectory(prefix='metaspici-reference-dryrun-') as td:
        tmp = Path(td)
        mags = tmp / 'mags'
        mags.mkdir()
        (mags / 'dryrun_host.fa').write_text('>dryrun_host_contig\nACGTACGTACGT\n')
        mge = tmp / 'dryrun_mge.fa'
        mge.write_text('>dryrun_mge_contig\nTTGCAATTGCAA\n')
        ref_out = tmp / 'reference_out'
        command = ref_command(args, ref_out, mags, mge, extra, dry_run=True)
        print('Shared reference-stage dry-run command:', ' '.join(command), flush=True)
        subprocess.run(command, check=True, env=env)


def validate_assembly_outputs(stage: Path, sample: str):
    mag_dir = stage / '02_binning' / f'{sample}_final_binset'
    mge_fasta = stage / '05_mge_host_linking' / f'metator_mge_{sample}' / 'mges_binned.fa'
    helper = ROOT / 'scripts' / 'assembly_utils.py'
    subprocess.run([sys.executable, str(helper), 'validate-inputs', str(mag_dir), str(mge_fasta)], check=True)
    return mag_dir, mge_fasta


def write_manifest(stage: Path, sample: str, variant: str = "short_hic"):
    mag_dir = stage / '02_binning' / f'{sample}_final_binset'
    mge_fasta = stage / '05_mge_host_linking' / f'metator_mge_{sample}' / 'mges_binned.fa'
    paths = [stage / '01_assembly' / f'{sample}_final.contigs.fa', mag_dir, mge_fasta,
             stage / '02_binning' / f'{sample}_final_shortbinset',
             stage / '02_binning' / f'{sample}_final_hicbinset']
    if variant in ('long_hic', 'hybrid_hic'):
        paths.append(stage / '02_binning' / f'{sample}_large_contigs')
    out = stage / 'run_info' / 'assembly_manifest.tsv'
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('w') as handle:
        handle.write('artifact\tpath\tstatus\tsize_bytes\n')
        for path in paths:
            if path.is_dir():
                files = [p for p in path.rglob('*') if p.is_file()]
                size = sum(p.stat().st_size for p in files)
                status = 'present' if files else 'empty'
            else:
                size = path.stat().st_size if path.is_file() else 0
                status = 'present' if path.is_file() and size else 'missing_or_empty'
            handle.write(f'{path.name}\t{path.relative_to(stage)}\t{status}\t{size}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--env-root', required=True)
    parser.add_argument('--checkm-db', required=True)
    parser.add_argument('--checkv-db', required=True)
    parser.add_argument('--genomad-db', required=True)
    parser.add_argument('--basalt-prefix', required=True)
    parser.add_argument('--cache-dir', required=True)
    parser.add_argument('--conda-executable', default=shutil.which('conda') or 'conda')
    parser.add_argument('--short-r1')
    parser.add_argument('--short-r2')
    parser.add_argument('--long-reads')
    parser.add_argument('--long-type', choices=LONG_TYPES)
    parser.add_argument('--opera-ms')
    parser.add_argument('--assembly-variant', choices=('short_hic','long_hic','hybrid_hic','hic_only'), default='short_hic')
    parser.add_argument('--hic-r1', required=True)
    parser.add_argument('--hic-r2', required=True)
    parser.add_argument('--enzyme', choices=['Sau3AI', 'MboI', 'DpnII'], required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--threads', type=int, default=32)
    parser.add_argument('--identity', type=float, default=99)
    parser.add_argument('--dry-run', action='store_true')
    args, ref_extra = parser.parse_known_args()
    args.reference_args = ref_extra
    try:
        check_inputs(args)
        if any(x in ('--checkm-db', '--checkv-db', '--genomad-db') or x.startswith(('--checkm-db=', '--checkv-db=', '--genomad-db=')) for x in ref_extra):
            raise ValueError('Database paths are configured with metaspici setup and cannot be overridden per run in assembly mode.')
        if any(x == '--re1' or x.startswith('--re1=') for x in ref_extra):
            raise ValueError('Use --enzyme only in assembly mode; do not also pass --re1.')
        out = Path(args.outdir).expanduser().resolve()
        out.mkdir(parents=True, exist_ok=True)
        stage, cfg_path, stage_cfg = prepare_config(args, out)
        env = conda_environment_variables(args.env_root)
        env['CHECKM_DATA_PATH'] = args.checkm_db
        env['CHECKM_DATA_DIR'] = args.checkm_db
        env['CHECKVDB'] = args.checkv_db
        run_snakemake(stage, cfg_path, args, env)
        if args.dry_run:
            reference_dry_run(args, out, ref_extra, env)
            print('Dry-run checks completed for the assembly DAG and shared reference DAG; no scientific jobs were executed.')
            return
        mag_dir, mge_fasta = validate_assembly_outputs(stage, args.sample)
        write_manifest(stage, args.sample, args.assembly_variant)
        command = ref_command(args, out, mag_dir, mge_fasta, ref_extra, dry_run=False)
        print('Starting shared reference workflow:', ' '.join(command), flush=True)
        subprocess.run(command, check=True, env=env)
        print(f'Assembly mode completed. Numbered assembly intermediates: {stage}; final reference results: {out / "results"}', flush=True)
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        parser.exit(1, f'metaspici assembly: {exc}\n')


if __name__ == '__main__':
    main()
