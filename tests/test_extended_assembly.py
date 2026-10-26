"""Regression checks for three new assembly variants (offline, no HPC executables)."""
import argparse
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from metaspici import cli
from metaspici.workflow import metaspici_assembly as runner
from metaspici.workflow.scripts import assembly_utils as util


class ExtendedAssemblyTest(unittest.TestCase):
    def _args(self, **kw):
        values = dict(short_r1=None, short_r2=None, long_reads=None, long_type=None,
                      opera_ms=None, hic_r1='/data/hic_R1.fq', hic_r2='/data/hic_R2.fq',
                      enzyme='Sau3AI', sample='X', outdir='/out', threads=8,
                      identity=99, dry_run=True)
        values.update(kw)
        return argparse.Namespace(**values)

    def test_cli_dispatches_all_four_input_patterns(self):
        patterns = [
            ({'short_r1':'/data/s1', 'short_r2':'/data/s2'}, 'short_hic'),
            ({'long_reads':'/data/long', 'long_type':'pacbioHiFi'}, 'long_hic'),
            ({'short_r1':'/data/s1', 'short_r2':'/data/s2',
              'long_reads':'/data/long', 'long_type':'nanopore',
              'opera_ms':'/some/OPERA-MS.pl'}, 'hybrid_hic'),
            ({}, 'hic_only'),
        ]
        config = dict(env_root='/envs', checkm_db='/db/checkm', checkv_db='/db/checkv',
                      genomad_db='/db/genomad', basalt_prefix='/envs/BASALT',
                      cache_dir='/cache', setup_mode='assembly')
        for kwargs, expected in patterns:
            with self.subTest(expected=expected):
                with patch.object(cli, 'validate_databases'), \
                     patch.object(cli, 'check_basalt'), \
                     patch.object(cli.shutil, 'which', return_value='/usr/bin/conda'), \
                     patch.object(cli.subprocess, 'run') as call:
                    cli.run_assembly(config, self._args(**kwargs), [])
                cmd = call.call_args.args[0]
                self.assertEqual(cmd[cmd.index('--assembly-variant') + 1], expected)
                self.assertEqual('--long-reads' in cmd, bool(kwargs.get('long_reads')))
                self.assertEqual('--short-r1' in cmd, bool(kwargs.get('short_r1')))
                self.assertEqual(cmd[cmd.index('--enzyme') + 1], 'Sau3AI')

    def test_incomplete_pairs_and_long_types_rejected(self):
        cfg = dict(env_root='/envs', checkm_db='/db/checkm', checkv_db='/db/checkv',
                   genomad_db='/db/genomad', basalt_prefix='/envs/BASALT',
                   cache_dir='/cache', setup_mode='assembly')
        with patch.object(cli, 'validate_databases'), patch.object(cli, 'check_basalt'):
            with self.assertRaisesRegex(ValueError, 'both --short-r1'):
                cli.run_assembly(cfg, self._args(short_r1='/r1'), [])
            with self.assertRaisesRegex(ValueError, 'requires --long-type'):
                cli.run_assembly(cfg, self._args(long_reads='/long'), [])
            with self.assertRaisesRegex(ValueError, 'needs --opera-ms'):
                cli.run_assembly(cfg, self._args(short_r1='/r1', short_r2='/r2',
                                               long_reads='/long', long_type='nanopore'), [])

    def test_split_long_contigs_and_merge_three_sources(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            long_in = base / 'contigs.fa'
            long_in.write_text('>big\n' + 'A' * 1000000 + '\n>small\nTTTT\n')
            singles = base / 'long'; self.assertEqual(util.split_long_contigs(long_in, singles), 1)
            self.assertEqual([p.name for p in singles.iterdir()], ['long_000001.fa'])
            hic = base / 'hic'; short = base / 'short'
            hic.mkdir(); short.mkdir()
            (hic / 'bin1.fa').write_text('>hic\nACGT\n')
            (short / 'bin1.fa').write_text('>short\nCCCC\n')
            out = base / 'merge'
            self.assertEqual(util.combine_bins_many(out, {'hic':hic, 'short':short, 'large':singles}), 3)
            self.assertEqual(sorted(p.name for p in out.iterdir()),
                             ['hic_000001.fa','large_000001.fa','short_000001.fa'])

    def test_truncater_real_fastq_file_naming(self):
        """The HiCUP output name used in Snakemake must match the real Perl tool."""
        import subprocess
        trunc = cli.WORKFLOW / 'third_party/hicup/hicup_truncater'
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp); out = base/'out'; out.mkdir()
            for n in (1,2):
                seq = 'AAAACCCCGATCGATCTTTTGGGG'
                (base/f'S_hic_clean_{n}.fastq').write_text(f'@r/{n}\n{seq}\n+\n'+ 'I'*len(seq)+'\n')
            result = subprocess.run(['perl',str(trunc),'--r','0','--re1','^GATC,Sau3AI',
                '--threads','1','--outdir',str(out),str(base/'S_hic_clean_1.fastq'),
                str(base/'S_hic_clean_2.fastq')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            for n in (1,2):
                self.assertTrue((out/f'S_hic_clean_{n}.trunc.fastq').is_file())


    def test_hybrid_uses_configured_opera_ms_without_per_run_override(self):
        config = dict(env_root='/envs', checkm_db='/db/checkm', checkv_db='/db/checkv',
                      genomad_db='/db/genomad', basalt_prefix='/envs/BASALT',
                      cache_dir='/cache', setup_mode='assembly',
                      opera_ms='/envs/opera_ms/OPERA-MS/OPERA-MS.pl')
        args = self._args(short_r1='/data/s1', short_r2='/data/s2',
                          long_reads='/data/long', long_type='nanopore', opera_ms=None)
        with patch.object(cli, 'validate_databases'), patch.object(cli, 'check_basalt'), \
             patch.object(cli.shutil, 'which', return_value='/usr/bin/conda'), \
             patch.object(cli.subprocess, 'run') as call:
            cli.run_assembly(config, args, [])
        cmd = call.call_args.args[0]
        self.assertEqual(cmd[cmd.index('--opera-ms') + 1], config['opera_ms'])

    def test_opera_ms_requires_complete_tree_not_just_pl_file(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'OPERA-MS.pl').write_text('#!/usr/bin/env perl\n')
            with self.assertRaisesRegex(ValueError, 'installation is incomplete'):
                cli._opera_ms_ready(root)

    def test_opera_ms_setup_build_and_dependency_check_commands(self):
        from subprocess import CompletedProcess
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / 'OPERA-MS'
            for rel in ['OPERA-MS.pl', 'Makefile', 'bin/test_time.pl', 'OPERA-LG', 'SIGMA',
                        'SHORT_READ_ANALYSIS', 'tools_opera_ms/install_perl_module.pl',
                        'tools_opera_ms/install_mummer3.23.sh']:
                p = root / rel; p.parent.mkdir(parents=True, exist_ok=True); p.write_text('test\n')
            with patch.object(cli.subprocess, 'run', side_effect=[
                    CompletedProcess(args=['make'], returncode=0),
                    CompletedProcess(args=['check-dependency'], returncode=0, stdout='', stderr='')
                ]) as call:
                result = cli.install_opera_ms('/fake/opera_ms', root, reuse=False)
            self.assertEqual(result, root / 'OPERA-MS.pl')
            self.assertTrue((root / '.metaspici_install_ok').is_file())
            self.assertIn('make', call.call_args_list[0].args[0])
            self.assertIn('check-dependency', call.call_args_list[1].args[0])

    def test_opera_ms_is_separated_from_quickmerge_environment(self):
        from pathlib import Path
        ref_engine = (cli.WORKFLOW / 'metaspici_ref.py').read_text()
        snakefile = (cli.WORKFLOW / 'Snakefile_assembly_extended').read_text()
        self.assertIn("'opera_ms': [['perl', '-v']", ref_engine)
        self.assertIn("'hybrid_tools': [['python', '-c', 'import shutil; assert shutil.which(\"merge_wrapper.py\")']]", ref_engine)
        self.assertIn("conda: 'opera_ms'", snakefile)
        self.assertIn("cfg['opera_ms_env'] = str(opera_prefix)", (cli.WORKFLOW.parent / 'cli.py').read_text())

    def test_opera_ms_dependency_check_uses_dedicated_prefix(self):
        from subprocess import CompletedProcess
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            install = base / 'envs' / 'opera_ms' / 'OPERA-MS'
            for rel in ['OPERA-MS.pl', 'Makefile', 'bin/test_time.pl', 'OPERA-LG', 'SIGMA',
                        'SHORT_READ_ANALYSIS', 'tools_opera_ms/install_perl_module.pl',
                        'tools_opera_ms/install_mummer3.23.sh']:
                target = install / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text('test\n')
            cfg = {'setup_mode': 'assembly', 'env_root': str(base / 'envs'),
                   'opera_ms_env': str(base / 'envs' / 'opera_ms'),
                   'opera_ms': str(install / 'OPERA-MS.pl')}
            with patch.object(cli, 'conda_run', return_value=CompletedProcess(args=[], returncode=0,
                                                                               stdout='', stderr='')) as run:
                cli.check_opera_ms(cfg)
            self.assertEqual(Path(run.call_args.args[0]), base / 'envs' / 'opera_ms')
            self.assertIn('check-dependency', run.call_args.args[1])

    def test_opera_ms_environment_declares_samtools_and_build_dependencies(self):
        import yaml
        env_root = cli.WORKFLOW / 'envs'
        opera = yaml.safe_load((env_root / 'opera_ms.yaml').read_text())
        hybrid = yaml.safe_load((env_root / 'hybrid_tools.yaml').read_text())
        self.assertIn('samtools=1.21', opera['dependencies'])
        self.assertIn('perl-app-cpanminus', opera['dependencies'])
        self.assertIn('mummer=3.23', opera['dependencies'])
        self.assertIn('quickmerge=0.3', hybrid['dependencies'])
        self.assertNotIn('samtools=1.21', hybrid['dependencies'])

    def test_original_short_workflow_unchanged(self):
        original = cli.WORKFLOW / 'Snakefile_assembly'
        extended = cli.WORKFLOW / 'Snakefile_assembly_extended'
        a = original.read_text(); e = extended.read_text()
        self.assertIn('rule short_assembly:', a)
        self.assertIn("if MODE == 'hic_only':", e)
        self.assertIn('rule hicup_for_assembly:', e)
        self.assertIn('rule opera_hybrid_assembly:', e)
        self.assertIn('merge_wrapper.py -l 1000', e)
        self.assertIn('rule extract_large_contigs:', e)
        self.assertIn('rule genomad_detection:', e)
        self.assertIn('rule metator_mge_link:', e)
        self.assertNotIn('hicstuff cutsite', e)


if __name__ == '__main__':
    unittest.main()
