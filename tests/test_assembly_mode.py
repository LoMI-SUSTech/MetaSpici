import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from metaspici import cli
from metaspici.workflow import metaspici_ref
from metaspici.workflow.scripts import assembly_utils


class AssemblyModeTest(unittest.TestCase):
    def test_assembly_mode_dispatches_and_preserves_common_options(self):
        argv = [
            'metaspici', 'run', '--mode', 'assembly',
            '--short-r1', 'short_R1.fastq.gz', '--short-r2', 'short_R2.fastq.gz',
            '--hic-r1', 'hic_R1.fastq.gz', '--hic-r2', 'hic_R2.fastq.gz',
            '--enzyme', 'Sau3AI', '--sample', 'sample1', '--threads', '16',
            '--identity', '97', '--outdir', 'out', '--min-depth', '0.16',
            '--fe-cutoff', '60', '--dry-run'
        ]
        with patch.object(sys, 'argv', argv), \
             patch.object(cli, 'read_config', return_value={'env_root': '/envs'}), \
             patch.object(cli.shutil, 'which', return_value='/usr/bin/conda'), \
             patch.object(cli, 'run_assembly') as run:
            cli.main()
        cfg, args, extra = run.call_args.args
        self.assertEqual(args.mode, 'assembly')
        self.assertEqual(args.identity, 97)
        self.assertEqual(args.threads, 16)
        self.assertIn('--min-depth', extra)
        self.assertEqual(extra[extra.index('--min-depth') + 1], '0.16')
        self.assertIn('--fe-cutoff', extra)
        self.assertTrue(args.dry_run)

    def test_profiles_keep_reference_environment_set(self):
        self.assertEqual(
            metaspici_ref.profile_envs('reference'),
            ('core', 'metaxa', 'hicup_runtime', 'bamm', 'genome_qc', 'viral')
        )
        self.assertTrue(set(metaspici_ref.profile_envs('reference')).issubset(
            set(metaspici_ref.profile_envs('assembly'))))
        self.assertIn('metator', metaspici_ref.profile_envs('assembly'))
        self.assertIn('genomad', metaspici_ref.profile_envs('assembly'))

    def test_combine_bins_prefixes_source_and_validates_fastas(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            hic = root / 'hic'; short = root / 'short'; combined = root / 'combined'
            hic.mkdir(); short.mkdir()
            (hic / 'bin1.fa').write_text('>c1\nACGT\n')
            (short / 'bin1.fa').write_text('>c2\nTTAA\n')
            count = assembly_utils.combine_bins(hic, short, combined)
            self.assertEqual(count, 2)
            self.assertTrue((combined / 'hic_bin1.fa').is_file())
            self.assertTrue((combined / 'short_bin1.fa').is_file())

    def test_resolve_basalt_binset_supports_legacy_and_release_names(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            legacy = root / 'Final_bestbinset'
            legacy.mkdir()
            (legacy / 'bin1.fa').write_text('>c1\nACGT\n')
            self.assertEqual(assembly_utils.resolve_basalt_binset(root), legacy)

    def test_validate_generated_inputs_rejects_empty_or_colliding_mge(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            mags = root / 'mags'; mags.mkdir()
            (mags / 'host.fa').write_text('>contig1\nACGT\n')
            mge = root / 'mge.fa'
            mge.write_text('>host\nTTAA\n')
            with self.assertRaisesRegex(ValueError, 'collide'):
                assembly_utils.validate_generated_inputs(mags, mge)

    def test_assembly_output_layout_is_stage_named_and_versioned(self):
        source = cli.WORKFLOW
        assembly = (source / 'Snakefile_assembly').read_text()
        runner = (source / 'metaspici_assembly.py').read_text()
        self.assertIn("stage = out / 'assembly_workflow'", runner)
        self.assertIn("stage / '01_assembly'", runner)
        self.assertIn("stage / '02_binning'", runner)
        self.assertIn("stage / '05_mge_host_linking'", runner)
        for folder in ('01_assembly/', '02_binning/', '03_mge_detection/', '04_hic_network/', '05_mge_host_linking/'):
            self.assertIn(folder, assembly)
        self.assertIn('legacy assembly_stage layout', runner)

    def test_old_v030_output_layout_is_rejected_to_prevent_mixing(self):
        from argparse import Namespace
        import tempfile
        from pathlib import Path
        with tempfile.TemporaryDirectory() as td:
            out = Path(td)
            old = out / 'assembly_stage'
            old.mkdir()
            (old / 'config.yaml').write_text('software_version: 0.3.0\n')
            args = Namespace()
            with self.assertRaisesRegex(ValueError, 'Use a new --outdir'):
                cli.WORKFLOW.joinpath('metaspici_assembly.py')
                from metaspici.workflow import metaspici_assembly
                metaspici_assembly.prepare_config(args, out)

    def test_assembly_workflow_uses_shared_truncater_not_hicstuff_cutsite(self):
        source = cli.WORKFLOW
        assembly = (source / 'Snakefile_assembly').read_text()
        reference = (source / 'Snakefile').read_text()
        self.assertIn('metator pipeline', assembly)
        self.assertIn('genomad end-to-end', assembly)
        self.assertIn('metator mge', assembly)
        self.assertNotIn('hicstuff cutsite', assembly)
        self.assertIn('rule hicup_truncater:', reference)
        self.assertIn('rule fastuniq:', reference)


if __name__ == '__main__':
    unittest.main()
