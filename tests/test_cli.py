import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from metaspici import cli

class InterfaceTest(unittest.TestCase):

    def test_checkm_database_can_be_bootstrapped_only_during_first_install(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            checkm = root / 'envs/genome_qc/checkm_data'
            checkv = root / 'checkv_db'; checkv.mkdir()
            genomad = root / 'genomad_db'; genomad.mkdir()
            cfg = {'checkm_db': str(checkm), 'checkv_db': str(checkv), 'genomad_db': str(genomad)}
            cli.validate_databases(cfg, 'assembly', allow_checkm_bootstrap=True)
            with self.assertRaisesRegex(ValueError, 'CheckM database missing/incomplete'):
                cli.validate_databases(cfg, 'assembly')

    def test_version_is_v11(self):
        from metaspici import __version__
        self.assertEqual(__version__, '1.1')

    def test_assembly_help_lists_assembly_inputs_without_reference_engine_help(self):
        argv = ['metaspici', 'run', '--mode', 'assembly', '--help']
        with patch.object(sys, 'argv', argv), patch.object(cli, 'read_config', return_value={}), \
             patch.object(cli, 'engine') as engine, contextlib.redirect_stdout(io.StringIO()) as output:
            cli.main()
        text = output.getvalue()
        self.assertIn('--short-r1', text)
        self.assertIn('--long-reads', text)
        self.assertIn('Hi-C only', text)
        self.assertNotIn('usage: metaspici_ref.py', text)
        engine.assert_not_called()
    def test_target_command_forwarding(self):
        cfg = {'env_root': '/existing/envs', 'checkm_db': '/db/checkm', 'checkv_db': '/db/checkv'}
        argv = ['metaspici', 'run', '--mode', 'reference', '--mag-dir', 'MAGs/',
                '--votu', 'viruses.fa', '--hic-r1', 'sample_R1.fastq.gz',
                '--hic-r2', 'sample_R2.fastq.gz', '--identity', '99',
                '--enzyme', 'Sau3AI', '--threads', '32', '--outdir', 'results/', '--dereplicate']
        with patch.object(sys, 'argv', argv), patch.object(cli, 'read_config', return_value=cfg), \
             patch.object(cli.shutil, 'which', return_value='/bin/conda'), patch.object(cli, 'engine') as run:
            cli.main()
        args = run.call_args.args[0]
        for flag, value in [('--sample', 'sample'), ('--hic-read1', 'sample_R1.fastq.gz'),
                            ('--output', 'results/'), ('--identity', '99'), ('--re1', '^GATC,Sau3AI'),
                            ('--env-root', '/existing/envs'), ('--checkv-db', '/db/checkv')]:
            self.assertEqual(args[args.index(flag)+1], value)
        self.assertIn('--dereplicate', args)


    def test_conda_environment_root_is_exported_for_named_environments(self):
        from metaspici.workflow.metaspici_ref import conda_environment_variables
        with tempfile.TemporaryDirectory() as td:
            env = conda_environment_variables(td, {
                'PATH': '/usr/bin',
                'CONDA_ENVS_PATH': '/existing/conda/envs:/other/envs',
            })
            self.assertEqual(env['CONDA_ENVS_PATH'], os.pathsep.join([
                str(Path(td).resolve()), '/existing/conda/envs', '/other/envs'
            ]))


    def test_dereplicate_is_optional_and_defaults_off(self):
        import yaml
        from metaspici.workflow import metaspici_ref
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            mag_dir = root/'mags'; mag_dir.mkdir()
            (mag_dir/'A.fa').write_text('>a\nACGT\n')
            mge = root/'mge.fa'; mge.write_text('>mge\nACGT\n')
            r1 = root/'r1.fastq'; r1.write_text('@r1\nACGT\n+\nIIII\n')
            r2 = root/'r2.fastq'; r2.write_text('@r2\nACGT\n+\nIIII\n')
            out = root/'out'
            argv = ['metaspici_ref.py', '--env-root', str(root/'envs'),
                    '--mag-dir', str(mag_dir), '--mge-fasta', str(mge),
                    '--hic-read1', str(r1), '--hic-read2', str(r2),
                    '--sample', 'test', '--output', str(out), '--dry-run']
            with patch.object(sys, 'argv', argv), \
                 patch.object(metaspici_ref.subprocess, 'run', return_value=type('R', (), {'returncode': 0})()), \
                 contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(SystemExit) as ex:
                    metaspici_ref.main()
            self.assertEqual(ex.exception.code, 0)
            cfg = yaml.safe_load((out/'run_info/config.yaml').read_text())
            self.assertFalse(cfg['dereplicate'])
            self.assertEqual(cfg['min_depth'], 0.15)
            self.assertEqual(cfg['mag_files'], [str((mag_dir/'A.fa').resolve())])

    def test_setup_reuse_does_not_install(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'config.json'
            argv = ['metaspici', 'setup', '--reuse-envs', '--env-root', td, '--checkv-db', td]
            with patch.object(sys, 'argv', argv), patch.object(cli, 'config_path', return_value=path), \
                 patch.object(cli.shutil, 'which', return_value='/bin/conda'), \
                 patch.object(cli, 'doctor') as doctor, patch.object(cli, 'engine') as engine:
                cli.main()
                engine.assert_not_called()
                doctor.assert_called_once()
            self.assertEqual(json.loads(path.read_text())['env_root'], td)

    def test_failed_setup_does_not_save(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)/'config.json'
            with patch.object(sys, 'argv', ['metaspici','setup','--reuse-envs','--env-root',td]), \
                 patch.object(cli, 'config_path', return_value=path), \
                 patch.object(cli.shutil, 'which', return_value='/bin/conda'), \
                 patch.object(cli, 'doctor', side_effect=ValueError('failed')), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit): cli.main()
            self.assertFalse(path.exists())

    def test_bundled_truncater_actual_function(self):
        def local_perl(prefix, args, **kwargs):
            return subprocess.run(args, **kwargs)
        with patch.object(cli, 'conda_run', side_effect=local_perl):
            cli.truncater_test(Path('/unused'))

if __name__ == '__main__':
    unittest.main()
