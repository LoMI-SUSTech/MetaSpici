import contextlib
import io
import tempfile
import unittest
from pathlib import Path
import sys
import subprocess
from unittest.mock import patch

from metaspici import cli
from metaspici.workflow import metaspici_ref

OLD = 'sp.check_output("pairtools", stderr=sp.STDOUT, shell=True)'
NEW = 'sp.check_output("pairtools --help", stderr=sp.STDOUT, shell=True)'


class MetaTORSetupTest(unittest.TestCase):
    def _source(self, td, contents=None):
        io = Path(td) / 'lib/python3.10/site-packages/metator/io.py'
        io.parent.mkdir(parents=True, exist_ok=True)
        io.write_text('import subprocess as sp\nprobe = ' + (contents or OLD) + '\n')
        return io

    def test_python_version_alias_is_deduplicated(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            real = self._source(td)
            alias_dir = root / 'lib/python3.1'
            alias_dir.parent.mkdir(parents=True, exist_ok=True)
            try:
                alias_dir.symlink_to(root / 'lib/python3.10', target_is_directory=True)
            except OSError:
                self.skipTest('directory symlinks unavailable on this platform')
            self.assertEqual(cli._metator_io_file(root), real.resolve())

    def test_two_distinct_metator_io_files_are_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            first = self._source(td)
            second = root / 'lib/python3.11/site-packages/metator/io.py'
            second.parent.mkdir(parents=True, exist_ok=True)
            second.write_text('different file')
            with self.assertRaisesRegex(ValueError, 'found 2'):
                cli._metator_io_file(root)

    def test_patch_backups_original_and_is_idempotent(self):
        with tempfile.TemporaryDirectory() as td:
            io = self._source(td)
            self.assertTrue(cli._fix_metator_pairtools_check(Path(td)))
            self.assertIn(NEW, io.read_text())
            backup = io.with_name('io.py.metaspici-original.bak')
            self.assertIn(OLD, backup.read_text())
            self.assertFalse(cli._fix_metator_pairtools_check(Path(td)))
            self.assertIn(OLD, backup.read_text())

    def test_unknown_metator_source_never_modified(self):
        with tempfile.TemporaryDirectory() as td:
            io = self._source(td, 'sp.check_output("pairtools --version")')
            with self.assertRaisesRegex(ValueError, 'Unexpected MetaTOR source'):
                cli._fix_metator_pairtools_check(Path(td))
            self.assertNotIn(NEW, io.read_text())
            self.assertFalse(io.with_name('io.py.metaspici-original.bak').exists())

    def test_reused_env_with_docopt_does_not_reinstall(self):
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / 'metator/conda-meta').mkdir(parents=True)
            (Path(td) / 'metator/conda-meta/history').touch()
            with patch.object(cli, 'conda_run') as run, \
                 patch.object(cli, '_fix_metator_pairtools_check') as fix, \
                 patch.object(cli, '_check_metator_dependencies') as validate, \
                 patch.object(cli.subprocess, 'run') as install:
                run.return_value.returncode = 0
                cli.prepare_metator_environment(td, reuse_envs=True)
                fix.assert_called_once()
                validate.assert_called_once()
                install.assert_not_called()

    def test_missing_docopt_repaired_when_not_reusing(self):
        with tempfile.TemporaryDirectory() as td:
            (Path(td) / 'metator/conda-meta').mkdir(parents=True)
            (Path(td) / 'metator/conda-meta/history').touch()
            with patch.object(cli, 'conda_run') as run, \
                 patch.object(cli, '_fix_metator_pairtools_check'), \
                 patch.object(cli, '_check_metator_dependencies'), \
                 patch.object(cli.subprocess, 'run') as install:
                run.return_value.returncode = 1
                cli.prepare_metator_environment(td, reuse_envs=False)
                self.assertIn('docopt=0.6.2', install.call_args.args[0])
                self.assertNotIn('pairtools', install.call_args.args[0])


    def test_first_assembly_setup_can_create_checkm_database(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td) / 'envs'
            db = root / 'genome_qc' / 'checkm_data'
            marker = db / 'genome_tree' / 'genome_tree.derep.txt'
            argv = ['metaspici_ref.py', '--env-root', str(root), '--create-envs-only',
                    '--profile', 'assembly', '--checkm-db', str(db)]

            def fake_conda_create(command, check, env):
                prefix = Path(command[command.index('--prefix') + 1])
                (prefix / 'conda-meta').mkdir(parents=True, exist_ok=True)
                (prefix / 'conda-meta' / 'history').touch()
                if prefix.name == 'genome_qc':
                    marker.parent.mkdir(parents=True, exist_ok=True)
                    marker.write_text('test database marker')
                    self.assertEqual(env['CHECKM_DATA_DIR'], str(db.resolve()))
                    self.assertEqual(env['CHECKM_DATA_PATH'], str(db.resolve()))
                return subprocess.CompletedProcess(command, 0)

            with patch.object(sys, 'argv', argv), \
                 patch.object(metaspici_ref.subprocess, 'run', side_effect=fake_conda_create) as run, \
                 contextlib.redirect_stdout(io.StringIO()):
                metaspici_ref.main()

            self.assertTrue(marker.is_file())
            self.assertEqual(run.call_count, len(metaspici_ref.profile_envs('assembly')))

    def test_wildcards_match_all_logs(self):
        source = (cli.WORKFLOW / 'Snakefile_assembly').read_text()
        import re
        rule_chunks = re.split(r'(?m)^rule ', source)[1:]
        self.assertEqual(len(rule_chunks), 14)
        for chunk in rule_chunks:
            if '{sample}' in chunk.split('    shell:', 1)[0].split('    run:', 1)[0]:
                outputs = re.search(r'(?m)^    output:(.*?)(?=^    [a-z_]+:|\Z)', chunk, re.S | re.M)
                if outputs and '{sample}' in outputs.group(1):
                    logs = re.search(r'(?m)^    log:(.*?)(?=^    [a-z_]+:|\Z)', chunk, re.S | re.M)
                    self.assertIsNotNone(logs)
                    for path in re.findall(r'logs/[^\'"]+\.log', logs.group(1)):
                        self.assertIn('{sample}', path, chunk[:40])


if __name__ == '__main__':
    unittest.main()
