#!/usr/bin/env python3
"""Reference-guided workflow engine; profile options also support assembly dependency setup."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import yaml

ROOT = Path(__file__).resolve().parent
REFERENCE_ENVS = ('core', 'metaxa', 'hicup_runtime', 'bamm', 'genome_qc', 'viral')
ASSEMBLY_ENVS = ('assembly', 'metator', 'genomad', 'long_assembly', 'hybrid_tools', 'opera_ms')
ENVS = REFERENCE_ENVS + ASSEMBLY_ENVS

def profile_envs(profile):
    if profile == 'reference':
        return REFERENCE_ENVS
    if profile == 'assembly':
        # Assembly mode also runs the shared reference workflow afterwards.
        return REFERENCE_ENVS + ASSEMBLY_ENVS
    if profile == 'all':
        return ENVS
    raise ValueError(f'Unknown setup profile: {profile}')


def conda_environment_variables(env_root, inherited=None):
    """Expose prefix-created named environments to Conda/Snakemake.

    Snakemake's ``conda: 'viral'`` syntax exports/activates by environment
    name. The environments in this package are created as prefixes under
    ``env_root``, so that root must be included in CONDA_ENVS_PATH.
    """
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
TOOLS = {
    'core': [['fastp', '--version'], ['fastuniq', '-h'], ['bwa'],
             ['samtools', '--version'], ['sambamba', '--version'],
             ['bedtools', '--version'], ['blastn', '-version'],
             ['barrnap', '--version'], ['tRNAscan-SE', '-h'], ['samblaster', '--version']],
    'metaxa': [['metaxa2', '--help']], 'hicup_runtime': [['perl', str(ROOT / 'third_party/hicup/hicup_truncater'), '--version']],
    'bamm': [['python', str(ROOT / 'scripts/bamm_filter_compat.py'), '--self-test']],
    'genome_qc': [['dRep', '-h'], ['checkm', '-h']],
    'viral': [['checkv', '--help'], ['anicalc', '--help'], ['aniclust', '--help']],
    'assembly': [['megahit', '--version'], ['seqmagick', '--version']],
    'metator': [['metator', '--help'], ['pairtools', '--help']],
    'genomad': [['genomad', '--help']],
    'long_assembly': [['flye', '--version'], ['python', '-c', 'import shutil; assert shutil.which("metaMDBG")'], ['seqkit', 'version']],
    'hybrid_tools': [['python', '-c', 'import shutil; assert shutil.which("merge_wrapper.py")']],
    'opera_ms': [['perl', '-v'], ['cpanm', '--version'], ['make', '--version'],
                 ['wget', '--version'], ['R', '--version'], ['java', '-version'],
                 ['samtools', '--version']],
}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--env-root', type=Path, default=ROOT / '.envs')
    p.add_argument('--create-envs-only', action='store_true')
    p.add_argument('--check-envs', action='store_true')
    p.add_argument('--profile', choices=['reference', 'assembly', 'all'], default='reference',
                   help='Environment profile; assembly includes reference and assembly-only tools')
    p.add_argument('--only-env', choices=ENVS, help='Install/check one environment')
    p.add_argument('--mag-dir', '--mag_dir', type=Path)
    p.add_argument('--votu', type=Path, help='Viral FASTA; required unless --mge-fasta is used')
    p.add_argument('--mge-fasta', type=Path, help='Already curated mixed MGE/plasmid FASTA; bypass viral-specific processing')
    p.add_argument('--checkv-db', type=Path)
    p.add_argument('--checkm-db', type=Path)
    p.add_argument('--genome-info', type=Path, help='Optional dRep CSV: genome,completeness,contamination; genome includes .fa')
    p.add_argument('--dereplicate', action='store_true',
                   help='Run dRep/CheckM QC and dereplication on all input .fa genomes; disabled by default')
    p.add_argument('-h1', '--hic-read1', type=Path)
    p.add_argument('-h2', '--hic-read2', type=Path)
    p.add_argument('-n', '--sample')
    p.add_argument('-j', '--threads', type=int, default=32)
    p.add_argument('-o', '--output', type=Path)
    p.add_argument('--identity', type=float, default=99)
    p.add_argument('--min-completeness', type=float, default=50)
    p.add_argument('--max-contamination', type=float, default=5)
    p.add_argument('--min-genome-length', type=int, default=50000)
    p.add_argument('--viral-coverage', type=float, default=85)
    p.add_argument('--re1', default='^GATC,Sau3AI')
    p.add_argument('--nofill', action='store_true')
    p.add_argument('--min-depth', type=float, default=0.15,
                   help='Fixed minimum Depth_Total (default: 0.15); not multiplied by read length')
    p.add_argument('--fe-cutoff', type=float, default=50,
                   help='Strict FE > cutoff for MGE-host final results (default: 50); all_results is unfiltered by FE')
    p.add_argument('--min-contacts', type=int, default=4, help='Inclusive minimum; original code used >3')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--target', help='Optional output filename or non-wildcard rule')
    a = p.parse_args()
    a.env_root = a.env_root.expanduser().resolve()
    selected = [a.only_env] if a.only_env else profile_envs(a.profile)
    if a.create_envs_only:
        a.env_root.mkdir(parents=True, exist_ok=True)
        install_env = {**os.environ, 'CONDA_CHANNEL_PRIORITY': 'strict'}
        if a.checkm_db:
            db = a.checkm_db.expanduser().resolve()
            marker = db / 'genome_tree/genome_tree.derep.txt'
            genome_qc_prefix = a.env_root / 'genome_qc'
            can_bootstrap_checkm = (
                'genome_qc' in selected
                and not (genome_qc_prefix / 'conda-meta/history').is_file()
            )
            if not marker.is_file() and not can_bootstrap_checkm:
                p.error('--checkm-db must contain genome_tree/genome_tree.derep.txt; '
                        'automatic first-time setup is only possible while creating a new genome_qc environment')
            # Bioconda post-link uses CHECKM_DATA_DIR; runtime uses CHECKM_DATA_PATH.
            # On a first install, CheckM may create/download the database at this path.
            install_env.update(CHECKM_DATA_DIR=str(db), CHECKM_DATA_PATH=str(db))
        for name in selected:
            prefix = a.env_root / name
            if (prefix / 'conda-meta/history').exists():
                missing = [cmd[0] for cmd in TOOLS[name] if not (prefix / 'bin' / cmd[0]).exists()]
                if missing:
                    p.error(f'Incomplete environment {prefix}: missing {missing}. Remove this prefix with conda env remove -p, then retry installation.')
                print(f'Existing environment retained: {prefix}', flush=True)
                continue
            subprocess.run(['conda', 'env', 'create', '--yes', '--prefix', str(prefix),
                            '--file', str(ROOT / 'envs' / f'{name}.yaml')], check=True,
                           env=install_env)
        if a.checkm_db:
            marker = a.checkm_db.expanduser().resolve() / 'genome_tree/genome_tree.derep.txt'
            if 'genome_qc' in selected and not marker.is_file():
                p.error(f'CheckM database was not found after environment setup: {marker}. '
                        'Check the CheckM installation log or provide an existing database path.')
        return
    if a.check_envs:
        failed = []
        env_for_names = conda_environment_variables(a.env_root)
        for name in selected:
            # Snakemake references these pre-created prefixes by name. Verify
            # Conda can resolve each name with the configured env-root first.
            exported = subprocess.run(
                ['conda', 'env', 'export', '--name', name],
                env=env_for_names, stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE, text=True)
            if exported.returncode != 0:
                detail = (exported.stderr or '').strip().splitlines()
                failed.append(f'{name}: Conda cannot resolve environment name under {a.env_root}'
                              + (f" ({detail[-1]})" if detail else ''))
            for command in TOOLS[name]:
                binary = a.env_root / name / 'bin' / command[0]
                if not binary.exists():
                    failed.append(f'{name}: missing {command[0]}')
                    continue
                result = subprocess.run(['conda', 'run', '--no-capture-output', '-p',
                                         str(a.env_root / name), *command])
                # BWA/FastUniq help can return 1 without indicating an installation failure.
                acceptable = (0, 1) if command[0] in ('bwa', 'fastuniq') else (0,)
                if result.returncode not in acceptable:
                    failed.append(f'{name}: {command} exited {result.returncode}')
        if failed:
            p.exit(1, '\n'.join(failed) + '\n')
        print('Executable and named-environment checks passed. Databases and real-data execution still require verification.')
        return
    for key in ('mag_dir', 'hic_read1', 'hic_read2', 'sample', 'output'):
        if getattr(a, key) is None:
            p.error(f'--{key.replace("_", "-")} is required')
    if bool(a.votu) == bool(a.mge_fasta):
        p.error('Supply exactly one of --votu or --mge-fasta')
    if a.votu and not a.checkv_db:
        p.error('--votu requires --checkv-db')
    if a.dereplicate and not a.checkm_db and not a.genome_info:
        p.error('--dereplicate requires --checkm-db, or existing CheckM results via --genome-info CSV')
    if not re.fullmatch(r'[A-Za-z0-9_-]+', a.sample):
        p.error('Sample must contain only letters, digits, underscore or hyphen')
    if not 95 <= a.identity <= 99:
        p.error('--identity must be between 95 and 99 (percent)')
    if a.threads < 1 or a.min_contacts < 1 or a.min_genome_length < 1:
        p.error('threads, min-contacts and min-genome-length must be positive')
    if not all(0 <= x <= 100 for x in (a.min_completeness, a.max_contamination, a.viral_coverage)):
        p.error('Quality and coverage percentages must lie in [0,100]')
    if a.min_depth < 0 or a.fe_cutoff < 0:
        p.error('Minimum depth and FE cutoff must be nonnegative')
    if not re.fullmatch(r'[ACGTN]*\^[ACGTN]+,[A-Za-z0-9_-]+', a.re1):
        p.error('--re1 must be a HiCUP enzyme definition, e.g. ^GATC,Sau3AI')
    for key in ('mag_dir', 'hic_read1', 'hic_read2', 'votu', 'mge_fasta',
                'checkm_db', 'checkv_db', 'genome_info'):
        if key in ('checkm_db', 'genome_info') and not a.dereplicate:
            continue
        value = getattr(a, key)
        if value:
            value = value.expanduser().resolve()
            if not value.exists():
                p.error(f'Not found: {value}')
            setattr(a, key, value)
    if not a.mag_dir.is_dir():
        p.error('--mag-dir must be a directory, one bacterial/archaeal genome per .fa file')
    mags = sorted(a.mag_dir.glob('*.fa'))
    if not mags:
        p.error('No *.fa genomes in --mag-dir; rename .fna/.fasta to .fa first')
    for f in mags:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', f.stem):
            p.error(f'Genome filename contains unsupported characters: {f.name}')
    out = a.output.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    cfg = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    cfg.update(mag_files=[str(f) for f in mags],
               root=str(ROOT), output=str(out),
               identity_fraction=a.identity / 100, sample=a.sample,
               dereplicate=bool(a.dereplicate))
    if (out / 'config.yaml').exists():
        p.error('This output directory uses the old layout; choose a new output directory')
    run_info = out / 'run_info'
    run_info.mkdir(exist_ok=True)
    cfg['output_layout_version'] = '0.2.2'
    cfg['software_version'] = '1.1'
    cf = run_info / 'config.yaml'
    # Refuse accidental reuse of an output directory with scientifically different settings.
    if cf.exists():
        old = yaml.safe_load(cf.read_text())
        ignored = {'dry_run', 'target', 'threads', 'check_envs', 'create_envs_only',
                   'only_env'}
        # Permit in-place upgrades from the previous validation releases. The new
        # dereplicate switch defaults to false; Snakemake will rebuild the host-genome
        # preparation output under the new rule logic.
        if old.get('software_version') in {'0.2.2', '0.2.3', '0.2.4', '0.2.5', '0.2.6'} and cfg.get('software_version') == '1.1':
            ignored.update({'software_version', 'drep_mag_files', 'preserve_genomes', 'dereplicate'})
        changed = [k for k in cfg if k not in ignored and old.get(k) != cfg[k]]
        if changed:
            p.error('Settings changed (' + ', '.join(changed) + '); use a new output directory')
    cf.write_text(yaml.safe_dump(cfg, sort_keys=False))
    command = ['snakemake', '-s', str(ROOT / 'Snakefile'), '--directory', str(out),
               '--configfile', str(cf), '--cores', str(a.threads), '--use-conda',
               '--printshellcmds', '--rerun-incomplete']
    if a.dry_run:
        command.append('--dry-run')
    if a.target:
        command.append(a.target)
    (run_info / 'invocation.json').write_text(json.dumps({'argv': sys.argv, 'command': command}, indent=2))
    # Ensure named Conda environments resolve to the already-created prefixes
    # under --env-root. Without this, Snakemake may treat an absolute prefix as
    # a name and call `conda env export --name /absolute/path`, which fails.
    run_env = conda_environment_variables(a.env_root)
    sys.exit(subprocess.run(command, env=run_env).returncode)

if __name__ == '__main__':
    main()
