"""Public command line for the MetaSpici reference and assembly workflows."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

from . import __version__

WORKFLOW = Path(__file__).resolve().parent / 'workflow'
TRUNCATER = WORKFLOW / 'third_party/hicup/hicup_truncater'
BASALT_REPOSITORY = 'https://github.com/PKU-EMBL/BASALT.git'
BASALT_REF = 'v1.1.0'  # legacy interface/output contract used by the original short_hic workflow
OPERA_MS_REPOSITORY = 'https://github.com/CSB5/OPERA-MS.git'
OPERA_MS_REF = 'OPERA-MS-0.9.0'


def config_path():
    return Path(os.environ.get('METASPICI_CONFIG', str(Path.home() / '.config/metaspici/config.json'))).expanduser()


def read_config():
    path = config_path()
    return json.loads(path.read_text()) if path.exists() else {}


def engine(args):
    subprocess.run([sys.executable, str(WORKFLOW / 'metaspici_ref.py'), *args], check=True)


def conda_run(prefix, args, **kwargs):
    return subprocess.run(['conda', 'run', '--no-capture-output', '-p', str(prefix), *args], **kwargs)


def truncater_test(root):
    with tempfile.TemporaryDirectory(prefix='metaspici-truncater-') as td:
        base = Path(td)
        out = base / 'out'
        out.mkdir()
        for mate in (1, 2):
            with (base / f'test_R{mate}.fastq').open('w') as f:
                for name, seq in [('junction', 'AAAACCCCGATCGATCTTTTGGGG'),
                                  ('control', 'AAAACCCCTTTTGGGG')]:
                    f.write(f'@{name}/{mate}\n{seq}\n+\n' + 'I' * len(seq) + '\n')
        conda_run(root / 'hicup_runtime', ['perl', str(TRUNCATER), '--re1', '^GATC,Sau3AI',
                  '--r', '0', '--threads', '1', '--outdir', str(out),
                  str(base / 'test_R1.fastq'), str(base / 'test_R2.fastq')], check=True)
        for mate in (1, 2):
            expected = (f'@junction/{mate}\nAAAACCCCGATC\n+\n' + 'I' * 12 + '\n' +
                        f'@control/{mate}\nAAAACCCCTTTTGGGG\n+\n' + 'I' * 16 + '\n')
            if (out / f'test_R{mate}.trunc.fastq').read_text() != expected:
                raise ValueError('Truncater functional check failed')
    print('PASS: HiCUP Truncater paired reads, sequences and qualities', flush=True)


def resolve_path(value):
    return str(Path(value).expanduser().resolve())


def validate_databases(cfg, mode='reference', allow_checkm_bootstrap=False):
    if cfg.get('checkm_db'):
        cm = Path(cfg['checkm_db'])
        if not (cm / 'genome_tree/genome_tree.derep.txt').is_file() and not allow_checkm_bootstrap:
            raise ValueError(f'CheckM database missing/incomplete: {cm}')
    if cfg.get('checkv_db') and not Path(cfg['checkv_db']).is_dir():
        raise ValueError(f"CheckV database directory not found: {cfg['checkv_db']}")
    if mode in ('assembly', 'all'):
        if not cfg.get('checkm_db'):
            raise ValueError('Assembly mode requires a CheckM database; set it with metaspici setup --checkm-db PATH')
        if not cfg.get('checkv_db'):
            raise ValueError('Assembly mode requires a CheckV database for MetaTOR; set it with metaspici setup --checkv-db PATH')
        if not cfg.get('genomad_db'):
            raise ValueError('Assembly mode requires a geNomad database; set it with metaspici setup --genomad-db PATH')
        for key, label in [('checkv_db', 'CheckV'), ('genomad_db', 'geNomad')]:
            db = Path(cfg[key])
            if not db.is_dir():
                raise ValueError(f'{label} database directory not found: {db}')


def _basalt_model_cache_ready():
    cache = Path.home() / '.cache'
    if (cache / 'BASALT.zip').is_file() and (cache / 'BASALT.zip').stat().st_size > 0:
        return True
    if cache.exists():
        return any(cache.rglob('*_ensemble.csv')) or any(cache.rglob('*.pt')) or any(cache.rglob('*.pth'))
    return False


def _basalt_source(cache_dir, allow_download=True):
    source = cache_dir / f'BASALT-{BASALT_REF}'
    required = ('basalt_env.yml', 'BASALT_script.zip', 'BASALT_models_download.py')
    if source.is_dir() and all((source / item).is_file() for item in required):
        return source
    if source.exists():
        raise ValueError(f'Incomplete BASALT source cache at {source}; move it aside and rerun setup')
    if not allow_download:
        raise ValueError(f'BASALT source not cached at {source}; rerun metaspici setup without --reuse-envs')
    if not shutil.which('git'):
        raise ValueError('git is required to install BASALT v1.1.0; install git or provide an existing --basalt-prefix')
    source.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(['git', 'clone', '--depth', '1', '--branch', BASALT_REF, BASALT_REPOSITORY, str(source)], check=True)
    missing = [item for item in required if not (source / item).is_file()]
    if missing:
        raise ValueError(f'BASALT {BASALT_REF} source is missing required files: {missing}')
    return source


def _install_basalt_scripts(source, prefix):
    from pathlib import PurePosixPath

    bin_dir = prefix / 'bin'
    bin_dir.mkdir(parents=True, exist_ok=True)
    archive_path = source / 'BASALT_script.zip'
    with zipfile.ZipFile(archive_path) as archive:
        members = [info for info in archive.infolist()
                   if info.filename and not info.filename.startswith('__MACOSX/')]
        top_levels = {PurePosixPath(info.filename).parts[0] for info in members
                      if PurePosixPath(info.filename).parts}
        strip_root = None
        if len(top_levels) == 1:
            candidate = next(iter(top_levels))
            if all(len(PurePosixPath(info.filename).parts) > 1 or info.is_dir() for info in members):
                strip_root = candidate
        for info in members:
            parts = PurePosixPath(info.filename).parts
            if not parts:
                continue
            if strip_root and parts[0] == strip_root:
                parts = parts[1:]
            if not parts:
                continue
            dest = bin_dir.joinpath(*parts)
            if info.is_dir():
                dest.mkdir(parents=True, exist_ok=True)
                continue
            dest.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(info) as src_handle, dest.open('wb') as dst_handle:
                shutil.copyfileobj(src_handle, dst_handle)
            mode = (info.external_attr >> 16) & 0o777
            if mode:
                dest.chmod(mode)
            if dest.name == 'BASALT' or dest.suffix in {'.py', '.pl'}:
                dest.chmod(dest.stat().st_mode | 0o111)
    # Some legacy script archives contain BASALT.py but no launcher.
    launcher = bin_dir / 'BASALT'
    if not launcher.exists() and (bin_dir / 'BASALT.py').is_file():
        launcher.write_text('#!/usr/bin/env bash\nexec "${CONDA_PREFIX}/bin/python" "${CONDA_PREFIX}/bin/BASALT.py" "$@"\n')
        launcher.chmod(0o755)
    if not launcher.is_file():
        raise ValueError(f'BASALT script archive did not provide a BASALT executable in {bin_dir}')


def _check_basalt_cli(prefix, env=None):
    probe = conda_run(prefix, ['BASALT', '--help'], capture_output=True, text=True, env=env)
    help_text = probe.stdout + probe.stderr
    if probe.returncode not in (0, 1) or ('BASALT' not in help_text) or ('--module' not in help_text) or ('--mode' not in help_text):
        raise ValueError(f'BASALT command check failed or lacks the legacy --mode/--module interface under {prefix}: {help_text[-1000:]}')


def install_basalt(prefix, cache_dir, checkm_db, reuse=False):
    """Install the legacy BASALT v1.1.0 CLI used by the original short_hic rules."""
    prefix = Path(prefix).expanduser().resolve()
    history = prefix / 'conda-meta/history'
    if reuse:
        if not history.is_file() or not (prefix / 'bin/BASALT').is_file():
            raise ValueError(f'Existing BASALT environment/executable not found at {prefix}; omit --reuse-envs to install it')
        if not _basalt_model_cache_ready():
            raise ValueError('BASALT model files were not found under ~/.cache; run setup without --reuse-envs to download them')
        _check_basalt_cli(prefix)
        return

    # Reuse a valid existing BASALT prefix and model cache even when the other
    # assembly environments still need installation.
    if history.is_file() and (prefix / 'bin/BASALT').is_file() and _basalt_model_cache_ready():
        _check_basalt_cli(prefix)
        return

    source = _basalt_source(Path(cache_dir).expanduser().resolve())
    env_vars = dict(os.environ)
    env_vars['CONDA_CHANNEL_PRIORITY'] = 'strict'
    env_vars['CHECKM_DATA_PATH'] = str(Path(checkm_db).expanduser().resolve())
    env_vars['CHECKM_DATA_DIR'] = str(Path(checkm_db).expanduser().resolve())
    if not history.is_file():
        prefix.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(['conda', 'env', 'create', '--yes', '--prefix', str(prefix),
                        '--file', str(source / 'basalt_env.yml')], check=True, env=env_vars)
    if not (prefix / 'bin/BASALT').is_file():
        _install_basalt_scripts(source, prefix)
    # The upstream v1.1.0 model downloader stores and loads weights from ~/.cache.
    # It skips network transfer when ~/.cache/BASALT.zip already exists.
    subprocess.run(['conda', 'run', '--no-capture-output', '-p', str(prefix),
                    'python', str(source / 'BASALT_models_download.py')],
                   check=True, cwd=str(source), env=env_vars)
    _check_basalt_cli(prefix, env=env_vars)


def check_basalt(cfg):
    if cfg.get('setup_mode') not in ('assembly', 'all'):
        return
    prefix = Path(cfg['basalt_prefix'])
    if not (prefix / 'bin/BASALT').is_file():
        raise ValueError(f'BASALT executable not found: {prefix}/bin/BASALT')
    if not _basalt_model_cache_ready():
        raise ValueError('BASALT model files were not found under ~/.cache. Run metaspici setup --mode assembly without --reuse-envs.')
    _check_basalt_cli(prefix)
    print(f'BASALT legacy --mode/--module interface check passed: {prefix}', flush=True)


def _opera_ms_ready(install_dir):
    """Check the complete upstream OPERA-MS tree, not just its Perl launcher."""
    root = Path(install_dir).expanduser().resolve()
    required = [
        root / 'OPERA-MS.pl', root / 'Makefile', root / 'bin' / 'test_time.pl',
        root / 'OPERA-LG', root / 'SIGMA', root / 'SHORT_READ_ANALYSIS',
        root / 'tools_opera_ms' / 'install_perl_module.pl',
        root / 'tools_opera_ms' / 'install_mummer3.23.sh',
    ]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise ValueError('OPERA-MS installation is incomplete; missing: ' + ', '.join(missing))
    return root


def _run_logged(command, cwd, log_path, env=None):
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open('w') as handle:
        result = subprocess.run(command, cwd=str(cwd), stdout=handle,
                                stderr=subprocess.STDOUT, env=env, text=True)
    if result.returncode:
        tail = log_path.read_text(errors='replace')[-4000:]
        raise ValueError(
            f'Command failed ({result.returncode}): {" ".join(map(str, command))}\n'
            f'Log: {log_path}\nLast log lines:\n{tail}'
        )
    return log_path


def install_opera_ms(prefix, install_dir, reuse=False):
    """Download/build OPERA-MS using its dedicated Conda environment.

    The environment prefix is env_root/opera_ms and the upstream checkout defaults
    to env_root/opera_ms/OPERA-MS. This keeps Perl, samtools, Java, R, compilers and
    OPERA-MS runtime checks isolated from core and QuickMerge environments.
    The optional reference database is not downloaded because MetaSpici uses
    --no-ref-clustering.
    """
    prefix = Path(prefix).expanduser().resolve()
    root = Path(install_dir).expanduser().resolve()
    if reuse:
        root = _opera_ms_ready(root)
        probe = conda_run(prefix, ['perl', str(root / 'OPERA-MS.pl'), 'check-dependency'],
                          capture_output=True, text=True)
        if probe.returncode:
            detail = (probe.stderr or probe.stdout or '').strip()[-3000:]
            raise ValueError('Existing OPERA-MS installation failed check-dependency. '
                             'Rerun setup without --reuse-envs to rebuild it. ' + detail)
        return root / 'OPERA-MS.pl'

    if not shutil.which('git'):
        raise ValueError('git is required for automatic OPERA-MS installation; install git and rerun setup.')
    if root.exists():
        try:
            root = _opera_ms_ready(root)
        except ValueError:
            # Recover only a partial checkout whose origin is the exact upstream repo.
            origin = subprocess.run(['git', '-C', str(root), 'remote', 'get-url', 'origin'],
                                    capture_output=True, text=True)
            if origin.returncode == 0 and origin.stdout.strip().rstrip('/') == OPERA_MS_REPOSITORY.rstrip('/'):
                shutil.rmtree(root)
            else:
                raise
    if not root.exists():
        root.parent.mkdir(parents=True, exist_ok=True)
        try:
            subprocess.run(['git', 'clone', '--depth', '1', '--branch', OPERA_MS_REF,
                            OPERA_MS_REPOSITORY, str(root)], check=True)
        except subprocess.CalledProcessError:
            # Avoid leaving an unusable partial upstream checkout behind.
            if root.exists():
                shutil.rmtree(root)
            raise

    root = _opera_ms_ready(root)
    marker = root / '.metaspici_install_ok'
    build_log = root / '.metaspici_make.log'
    dep_log = root / '.metaspici_check_dependency.log'
    env = dict(os.environ)
    env['CONDA_CHANNEL_PRIORITY'] = 'strict'
    if marker.is_file() and f'ref={OPERA_MS_REF}' in marker.read_text(errors='replace'):
        prior_probe = subprocess.run(['conda', 'run', '--no-capture-output', '-p', str(prefix),
                                      'perl', str(root / 'OPERA-MS.pl'), 'check-dependency'],
                                     cwd=str(root), capture_output=True, text=True, env=env)
        if prior_probe.returncode == 0:
            print(f'OPERA-MS already installed and dependency check passed: {root / "OPERA-MS.pl"}', flush=True)
            return root / 'OPERA-MS.pl'
    # The upstream Makefile builds OPERA-LG/SIGMA, installs bundled tools and
    # invokes tools_opera_ms/install_perl_module.pl. The dedicated opera_ms env
    # supplies Perl, cpanm, samtools, compiler, make, R and Java.
    _run_logged(['conda', 'run', '--no-capture-output', '-p', str(prefix), 'make'],
                root, build_log, env=env)
    probe = subprocess.run(['conda', 'run', '--no-capture-output', '-p', str(prefix),
                            'perl', str(root / 'OPERA-MS.pl'), 'check-dependency'],
                           cwd=str(root), capture_output=True, text=True, env=env)
    dep_log.write_text((probe.stdout or '') + (probe.stderr or ''))
    if probe.returncode:
        detail = dep_log.read_text(errors='replace')[-4000:]
        raise ValueError('OPERA-MS build completed but check-dependency failed. '
                         f'Log: {dep_log}\nLast log lines:\n{detail}')
    marker.write_text(f'repository={OPERA_MS_REPOSITORY}\nref={OPERA_MS_REF}\n')
    print(f'OPERA-MS installed and dependency check passed: {root / "OPERA-MS.pl"}', flush=True)
    print(f'OPERA-MS build log: {build_log}', flush=True)
    print(f'OPERA-MS dependency log: {dep_log}', flush=True)
    return root / 'OPERA-MS.pl'


def check_opera_ms(cfg):
    if cfg.get('setup_mode') not in ('assembly', 'all'):
        return
    value = cfg.get('opera_ms')
    if not value:
        raise ValueError('OPERA-MS is not configured. Rerun metaspici setup --mode assembly without --reuse-envs.')
    root = Path(value).expanduser().resolve()
    _opera_ms_ready(root.parent)
    prefix = Path(cfg.get('opera_ms_env', Path(cfg['env_root']) / 'opera_ms'))
    probe = conda_run(prefix, ['perl', str(root), 'check-dependency'], capture_output=True, text=True)
    if probe.returncode:
        detail = (probe.stderr or probe.stdout or '').strip()[-3000:]
        raise ValueError(f'OPERA-MS dependency check failed for {root}: {detail}')
    print(f'OPERA-MS dependency check passed: {root}', flush=True)


def doctor(cfg):
    root = Path(cfg['env_root'])
    profile = cfg.get('setup_mode', 'reference')
    engine(['--check-envs', '--profile', profile, '--env-root', str(root)])
    truncater_test(root)
    validate_databases(cfg, profile)
    check_basalt(cfg)
    if profile in ('assembly', 'all'):
        _check_metator_dependencies(root / 'metator')
        check_opera_ms(cfg)
    print('Environment checks passed. Database completeness and real-data execution still require validation.')


def _metator_io_file(prefix):
    """Find MetaTOR's io.py, deduplicating Python-version directory aliases.

    Some Conda environments expose ``lib/python3.1`` as an alias to
    ``lib/python3.10``. A glob using ``python*`` then finds the same file via
    two lexical paths. Resolve matches and deduplicate by canonical path before
    deciding whether the installation is ambiguous.
    """
    prefix = Path(prefix).resolve()
    candidates = sorted((prefix / 'lib').glob('python*/site-packages/metator/io.py'))
    matches = {}
    for candidate in candidates:
        try:
            resolved = candidate.resolve(strict=True)
            resolved.relative_to(prefix)
        except (OSError, ValueError):
            continue
        if resolved.is_file():
            matches[str(resolved)] = resolved
    unique_matches = [matches[key] for key in sorted(matches)]
    if len(unique_matches) != 1:
        found = '\n'.join(f'  - {p}' for p in unique_matches) or '  (none)'
        raise ValueError(
            f'Expected one unique MetaTOR io.py under {prefix}, found {len(unique_matches)}. '
            f'Candidates after resolving aliases:\n{found}; unable to apply the compatibility fix safely'
        )
    return unique_matches[0]


def _fix_metator_pairtools_check(prefix):
    """Apply the validated MetaTOR 1.3.7 CLI check fix once, with backup.

    MetaTOR 1.3.7 invokes bare `pairtools` to detect installation. Some
    pairtools versions return status 2 when run without a subcommand, leading
    to a false ImportError. This changes ONLY the dependency probe, not sorting.
    """
    import ast
    source = _metator_io_file(prefix)
    original = source.read_text()
    old = 'sp.check_output("pairtools", stderr=sp.STDOUT, shell=True)'
    new = 'sp.check_output("pairtools --help", stderr=sp.STDOUT, shell=True)'
    if original.count(new) == 1 and old not in original:
        print(f'MetaTOR pairtools check already compatible: {source}', flush=True)
        return False
    if original.count(old) != 1 or new in original:
        raise ValueError(f'Unexpected MetaTOR source in {source}; cannot safely apply '
                         'pairtools compatibility fix. No file was modified.')
    updated = original.replace(old, new, 1)
    ast.parse(updated, filename=str(source))
    backup = source.with_name(source.name + '.metaspici-original.bak')
    if not backup.exists():
        backup.write_text(original)
    # Write using a same-directory temp file, then replace atomically.
    with tempfile.NamedTemporaryFile(mode='w', dir=source.parent,
                                     prefix=source.name + '.', suffix='.tmp',
                                     delete=False) as handle:
        handle.write(updated)
        temp = Path(handle.name)
    try:
        temp.chmod(source.stat().st_mode)
        temp.replace(source)
    finally:
        if temp.exists():
            temp.unlink()
    print(f'MetaTOR pairtools CLI compatibility fix applied: {source}', flush=True)
    return True


def _check_metator_dependencies(prefix):
    """Check executable behavior, not just 'metator --help'."""
    probes = [
        (['python', '-c',
          'import docopt; assert getattr(docopt, "__version__", None) == "0.6.2"'],
         'docopt=0.6.2'),
        (['pairtools', '--help'], 'pairtools --help'),
    ]
    for command, label in probes:
        result = conda_run(prefix, command, capture_output=True, text=True)
        if result.returncode != 0:
            detail = (result.stderr or result.stdout or '').strip()[-1000:]
            raise ValueError(f'MetaTOR dependency check failed ({label}) in {prefix}: {detail}')
    print(f'MetaTOR docopt and pairtools checks passed: {prefix}', flush=True)


def prepare_metator_environment(root, reuse_envs=False):
    """Setup-only adjustment for assembly/all profiles; no reference changes."""
    prefix = Path(root) / 'metator'
    if not (prefix / 'conda-meta/history').is_file():
        raise ValueError(f'MetaTOR environment not installed: {prefix}')
    docopt_probe = conda_run(prefix, ['python', '-c',
                                      'import docopt; assert getattr(docopt, "__version__", None) == "0.6.2"'],
                             capture_output=True, text=True)
    if docopt_probe.returncode:
        if reuse_envs:
            raise ValueError('The retained MetaTOR environment lacks docopt=0.6.2. '
                             'Rerun metaspici setup --mode assembly without --reuse-envs '
                             'to repair the dependency automatically.')
        # Existing environments are normally skipped by --create-envs-only.
        # Repair docopt in-place without reinstalling pairtools or MetaTOR.
        subprocess.run(['conda', 'install', '-y', '-p', str(prefix),
                        '--freeze-installed', '--override-channels',
                        '-c', 'conda-forge', '-c', 'bioconda', 'docopt=0.6.2'], check=True)
    _fix_metator_pairtools_check(prefix)
    _check_metator_dependencies(prefix)


def _save_config(cfg):
    path = config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as handle:
        json.dump(cfg, handle, indent=2)
        tmp = handle.name
    os.replace(tmp, path)
    print(f'Setup complete. Configuration: {path}')


def _setup(args, cfg):
    root = Path(args.env_root or cfg.get('env_root', Path.home() / '.local/share/metaspici/envs')).expanduser().resolve()
    mode = args.mode
    cfg['env_root'] = str(root)
    cfg['setup_mode'] = mode
    if args.checkm_db:
        cfg['checkm_db'] = resolve_path(args.checkm_db)
    if args.checkv_db:
        cfg['checkv_db'] = resolve_path(args.checkv_db)
    if args.genomad_db:
        cfg['genomad_db'] = resolve_path(args.genomad_db)
    if args.cache_dir:
        cfg['cache_dir'] = resolve_path(args.cache_dir)
    else:
        cfg.setdefault('cache_dir', str(Path.home() / '.cache' / 'metaspici'))
    if args.basalt_prefix:
        cfg['basalt_prefix'] = resolve_path(args.basalt_prefix)
    else:
        cfg.setdefault('basalt_prefix', str(root / 'basalt'))
    opera_prefix = root / 'opera_ms'
    if args.opera_ms_dir:
        # Retain an advanced override for users with an existing upstream checkout;
        # the Conda runtime environment is still managed under env_root/opera_ms.
        cfg['opera_ms_dir'] = resolve_path(args.opera_ms_dir)
    else:
        cfg['opera_ms_dir'] = str(opera_prefix / 'OPERA-MS')
    cfg['opera_ms_env'] = str(opera_prefix)
    if mode in ('assembly', 'all'):
        checkm_db = Path(cfg['checkm_db']).expanduser().resolve() if cfg.get('checkm_db') else None
        checkm_marker = (checkm_db / 'genome_tree/genome_tree.derep.txt') if checkm_db else None
        genome_qc_history = root / 'genome_qc' / 'conda-meta' / 'history'
        allow_checkm_bootstrap = bool(
            not args.reuse_envs
            and checkm_db is not None
            and checkm_marker is not None
            and not checkm_marker.is_file()
            and not genome_qc_history.is_file()
        )
        validate_databases(cfg, mode, allow_checkm_bootstrap=allow_checkm_bootstrap)
    if not args.reuse_envs:
        install_args = ['--create-envs-only', '--profile', mode, '--env-root', str(root)]
        if cfg.get('checkm_db'):
            install_args += ['--checkm-db', cfg['checkm_db']]
        engine(install_args)
        if mode in ('assembly', 'all'):
            # A first-time CheckM Conda install may create/download checkm_data via
            # CHECKM_DATA_DIR/CHECKM_DATA_PATH; verify the expected marker afterwards.
            validate_databases(cfg, mode)
        # Keep the existing targeted setuptools repair for the legacy CheckM env.
        probe = conda_run(root / 'genome_qc', ['python', '-c', 'import pkg_resources'], capture_output=True, text=True)
        if probe.returncode and "No module named 'pkg_resources'" in probe.stderr:
            subprocess.run(['conda', 'install', '-p', str(root / 'genome_qc'),
                            '--override-channels', '-c', 'conda-forge', '-c', 'bioconda',
                            'setuptools=80.9.0', '-y'], check=True)
        if mode in ('assembly', 'all'):
            install_basalt(cfg['basalt_prefix'], cfg['cache_dir'], cfg['checkm_db'], reuse=False)
    if mode in ('assembly', 'all'):
        prepare_metator_environment(root, reuse_envs=args.reuse_envs)
        cfg['opera_ms'] = str(install_opera_ms(opera_prefix, cfg['opera_ms_dir'], reuse=args.reuse_envs))
    doctor(cfg)
    _save_config(cfg)


def run_assembly(cfg, args, extra):
    if not cfg.get('env_root'):
        raise ValueError('Run metaspici setup --mode assembly first.')
    if cfg.get('setup_mode') not in ('assembly', 'all'):
        raise ValueError('Assembly dependencies are not configured. Run metaspici setup --mode assembly with the required database paths.')
    validate_databases(cfg, 'assembly')
    check_basalt(cfg)
    if not (args.hic_r1 and args.hic_r2):
        raise ValueError('Assembly mode requires paired --hic-r1 and --hic-r2.')
    if bool(args.short_r1) != bool(args.short_r2):
        raise ValueError('Provide both --short-r1 and --short-r2, or neither.')
    if args.long_reads and not args.long_type:
        raise ValueError('Long-read input requires --long-type nanopore|pacbio|nanoporeR10|pacbioHiFi.')
    if args.long_type and not args.long_reads:
        raise ValueError('--long-type requires --long-reads.')
    if args.short_r1 and args.long_reads:
        variant = 'hybrid_hic'
    elif args.short_r1:
        variant = 'short_hic'
    elif args.long_reads:
        variant = 'long_hic'
    else:
        variant = 'hic_only'
    if variant == 'hybrid_hic':
        args.opera_ms = args.opera_ms or cfg.get('opera_ms')
        if not args.opera_ms:
            raise ValueError('Hybrid assembly needs --opera-ms or a configured OPERA-MS installation. Run metaspici setup --mode assembly without --reuse-envs to install it automatically.')
    if not args.enzyme:
        raise ValueError('Assembly mode requires --enzyme (Sau3AI, MboI, or DpnII).')
    if any(x == '--re1' or x.startswith('--re1=') for x in extra):
        raise ValueError('Assembly mode uses --enzyme for both MetaTOR and the shared reference workflow; do not also pass --re1.')
    if any(x in ('--checkm-db', '--checkv-db', '--genomad-db') or x.startswith(('--checkm-db=', '--checkv-db=', '--genomad-db=')) for x in extra):
        raise ValueError('Configure CheckM, CheckV, and geNomad database paths with metaspici setup, not per-run overrides.')
    if any(x in ('--mag-dir', '--votu', '--mge-fasta', '--output', '--hic-read1', '--hic-read2') or
           x.startswith(('--mag-dir=', '--votu=', '--mge-fasta=', '--output=', '--hic-read1=', '--hic-read2=')) for x in extra):
        raise ValueError('Assembly mode creates its own MAG and MGE inputs; do not pass reference-mode input path options.')
    command = [sys.executable, str(WORKFLOW / 'metaspici_assembly.py'),
               '--env-root', cfg['env_root'], '--checkm-db', cfg['checkm_db'],
               '--checkv-db', cfg['checkv_db'], '--genomad-db', cfg['genomad_db'],
               '--basalt-prefix', cfg['basalt_prefix'], '--cache-dir', cfg['cache_dir'],
               '--conda-executable', shutil.which('conda') or 'conda',
               '--assembly-variant', variant,
               '--hic-r1', args.hic_r1, '--hic-r2', args.hic_r2,
               '--enzyme', args.enzyme, '--sample', args.sample,
               '--outdir', args.outdir or '.', '--threads', str(args.threads or 32),
               '--identity', f"{args.identity if args.identity is not None else 99:g}"]
    if args.short_r1:
        command += ['--short-r1', args.short_r1, '--short-r2', args.short_r2]
    if args.long_reads:
        command += ['--long-reads', args.long_reads, '--long-type', args.long_type]
    if args.opera_ms:
        command += ['--opera-ms', args.opera_ms]
    if args.dry_run:
        command.append('--dry-run')
    command.extend(extra)
    subprocess.run(command, check=True)


def main():
    p = argparse.ArgumentParser(prog='metaspici')
    p.add_argument('--version', action='version', version=__version__)
    commands = p.add_subparsers(dest='command', required=True)
    setup = commands.add_parser('setup', help='Create/reuse environments, validate and save database paths')
    setup.add_argument('--mode', choices=['reference', 'assembly', 'all'], default='reference',
                       help='reference installs the existing profile; assembly also installs assembly tools and builds OPERA-MS for Hybrid mode')
    setup.add_argument('--env-root', type=Path)
    setup.add_argument('--checkm-db', type=Path)
    setup.add_argument('--checkv-db', type=Path)
    setup.add_argument('--genomad-db', type=Path)
    setup.add_argument('--basalt-prefix', type=Path, help='Existing or target prefix for the pinned legacy BASALT v1.1.0 environment')
    setup.add_argument('--opera-ms-dir', type=Path, help='Advanced override for an upstream OPERA-MS checkout; default: <env-root>/opera_ms/OPERA-MS')
    setup.add_argument('--cache-dir', type=Path, help='Source cache for BASALT v1.1.0 (default: ~/.cache/metaspici)')
    setup.add_argument('--reuse-envs', action='store_true', help='Do not install packages; validate existing environments')
    commands.add_parser('doctor', help='Check configured environments, databases and Truncater function')
    commands.add_parser('config', help='Show saved paths')
    run = commands.add_parser('run', help='Run reference-guided or assembly-based MGE-host inference', add_help=False)
    run.add_argument('--mode', choices=['reference', 'assembly'], default='reference')
    run.add_argument('--enzyme', choices=['Sau3AI', 'MboI', 'DpnII'])
    run.add_argument('--hic-r1')
    run.add_argument('--hic-r2')
    run.add_argument('--short-r1', help='Paired-end short-read R1 (assembly mode only)')
    run.add_argument('--short-r2', help='Paired-end short-read R2 (assembly mode only)')
    run.add_argument('--long-reads', help='ONT/PacBio FASTQ input (assembly mode only)')
    run.add_argument('--long-type', choices=['nanopore','pacbio','nanoporeR10','pacbioHiFi'])
    run.add_argument('--opera-ms', help='Optional override for the automatically installed OPERA-MS.pl (hybrid only)')
    run.add_argument('--outdir')
    run.add_argument('--sample')
    run.add_argument('-j', '--threads', type=int)
    run.add_argument('--identity', type=float)
    run.add_argument('--dry-run', action='store_true')
    run.add_argument('-h', '--help', action='store_true')
    args, extra = p.parse_known_args()
    if args.command not in ('run',) and extra:
        p.error('Unrecognized arguments: ' + ' '.join(extra))
    try:
        cfg = read_config()
        if args.command == 'config':
            print(json.dumps(cfg, indent=2)); return
        if getattr(args, 'help', False):
            run.print_help()
            if args.mode == 'reference':
                print('\nReference-mode options include --mag-dir, --votu/--mge-fasta, --min-depth, --fe-cutoff and --min-contacts.')
                engine(['--help'])
            else:
                print('\nAssembly input variants: Short reads + Hi-C; Long reads + Hi-C; Hybrid + Hi-C; Hi-C only.')
                print('Long-read routes accept --long-type nanopore, pacbio, nanoporeR10, or pacbioHiFi.')
            return
        if not shutil.which('conda'):
            raise ValueError('conda is not on PATH. Activate the launcher environment first.')
        if args.command == 'setup':
            _setup(args, cfg)
            return
        if args.command == 'doctor':
            if not cfg.get('env_root'):
                raise ValueError('Run metaspici setup first.')
            doctor(cfg)
            return
        if not cfg.get('env_root'):
            raise ValueError('Run metaspici setup first.')
        if args.enzyme and any(x == '--re1' or x.startswith('--re1=') for x in extra) and args.mode == 'reference':
            raise ValueError('Use either --enzyme or --re1, not both.')
        if args.mode == 'assembly':
            if not args.sample:
                sample_source = Path(args.hic_r1).name if args.hic_r1 else 'sample'
                sample_source = re.sub(r'\.(fastq|fq)(\.gz)?$', '', sample_source)
                sample_source = re.sub(r'([_.-]R?1)([_.-]001)?$', '', sample_source)
                args.sample = re.sub(r'[^A-Za-z0-9_-]', '_', sample_source) or 'sample'
            if not re.fullmatch(r'[A-Za-z0-9_-]+', args.sample):
                raise ValueError('Sample must contain only letters, digits, underscore or hyphen')
            if args.identity is not None and not 95 <= args.identity <= 99:
                raise ValueError('--identity must be between 95 and 99 (percent)')
            run_assembly(cfg, args, extra)
            return
        # Reference mode preserves the established reference workflow arguments and defaults.
        if args.short_r1 or args.short_r2 or args.long_reads or args.long_type or args.opera_ms:
            raise ValueError('Assembly read/tool options are only used with --mode assembly.')
        if args.identity is not None and not 95 <= args.identity <= 99:
            raise ValueError('--identity must be between 95 and 99 (percent)')
        forwarded = ['--env-root', cfg['env_root']]
        for key, flag in [('checkm_db', '--checkm-db'), ('checkv_db', '--checkv-db')]:
            if cfg.get(key) and not any(x == flag or x.startswith(flag+'=') for x in extra):
                forwarded += [flag, cfg[key]]
        for value, flag in [(args.hic_r1, '--hic-read1'), (args.hic_r2, '--hic-read2'), (args.outdir, '--output')]:
            if value:
                forwarded += [flag, value]
        if args.threads is not None:
            forwarded += ['--threads', str(args.threads)]
        if args.identity is not None:
            forwarded += ['--identity', f'{args.identity:g}']
        if args.dry_run:
            forwarded.append('--dry-run')
        sample = args.sample
        if not sample:
            sample = Path(args.hic_r1).name if args.hic_r1 else 'sample'
            sample = re.sub(r'\.(fastq|fq)(\.gz)?$', '', sample)
            sample = re.sub(r'([_.-]R?1)([_.-]001)?$', '', sample)
            sample = re.sub(r'[^A-Za-z0-9_-]', '_', sample) or 'sample'
        forwarded += ['--sample', sample]
        if args.enzyme:
            forwarded += ['--re1', f'^GATC,{args.enzyme}']
        engine(forwarded + extra)
    except (ValueError, OSError, subprocess.CalledProcessError, zipfile.BadZipFile) as exc:
        p.exit(1, f'metaspici: {exc}\n')


if __name__ == '__main__':
    main()
