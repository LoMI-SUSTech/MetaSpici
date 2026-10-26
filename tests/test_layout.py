from pathlib import Path
import tempfile
import unittest
import importlib.util
from metaspici.cli import WORKFLOW

spec = importlib.util.spec_from_file_location('organize', WORKFLOW/'scripts/organize_reports.py')
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
spec2 = importlib.util.spec_from_file_location('refutils', WORKFLOW/'scripts/ref_utils.py')
utils = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(utils)

class LayoutTest(unittest.TestCase):

    def test_snakefile_uses_named_precreated_conda_envs(self):
        snakefile = (WORKFLOW / 'Snakefile').read_text()
        self.assertIn("conda: 'viral'", snakefile)
        self.assertIn("CORE = 'core'", snakefile)
        self.assertNotIn("str(ENV / 'viral')", snakefile)


    def test_drep_is_optional_and_all_inputs_are_preserved_by_default(self):
        snakefile = (WORKFLOW / 'Snakefile').read_text()
        self.assertIn("dereplicate='true' if config.get('dereplicate', False) else 'false'", snakefile)
        self.assertIn("conda: 'genome_qc' if config.get('dereplicate', False) else CORE", snakefile)
        self.assertIn('if [ "{params.dereplicate}" = "true" ]; then', snakefile)
        self.assertIn('cp -- {input.genomes:q} {output:q}/dereplicated_genomes/', snakefile)
        self.assertIn("merge_genomes {input:q}/dereplicated_genomes", snakefile)
        self.assertNotIn("'yeast' in p.stem", (WORKFLOW / 'metaspici_ref.py').read_text())

        # Exercise the normal merge path with multiple arbitrary filenames. The
        # default workflow must not depend on a filename-based genome-type guess.
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            derep = root/'drep'/'dereplicated_genomes'
            derep.mkdir(parents=True)
            (derep/'A.fa').write_text('>a_contig\nACGT\n')
            (derep/'Eukaryote_01.fa').write_text('>euk_contig\nAACCGG\n')
            utils.merge_genomes(derep, str(root/'host_genomes.fasta'), str(root/'offsets.tsv'))
            merged = (root/'host_genomes.fasta').read_text()
            self.assertIn('>A\nACGT', merged)
            self.assertIn('>Eukaryote_01\nAACCGG', merged)

    def test_report_collection_and_manifest(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            files={'intermediates/reference_qc/drep/data_tables/Chdb.csv':'genome,completeness\nhost.fa,99\n',
                   'intermediates/reference_qc/checkv/quality_summary.tsv':'contig_id\nvirus\n',
                   'intermediates/reads/truncated/hicup_summary.txt':'summary\n',
                   'intermediates/reads/truncated/reads.fastq':'@x\nAC\n+\nII\n',
                   'results/sample.mge_host.filtered.tsv':'source\ttarget\nvirus\thost\n',
                   '.snakemake/log/internal.log':'internal'}
            for rel,text in files.items():
                p=root/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(text)
            module.organize(root,'sample')
            self.assertEqual((root/'qc/reference/tool_reports/drep/Chdb.csv').read_text(),files['intermediates/reference_qc/drep/data_tables/Chdb.csv'])
            self.assertTrue((root/'qc/reads/hicup/hicup_summary.txt').is_file())
            self.assertFalse((root/'qc/reads/hicup/reads.fastq').exists())
            inventory=(root/'run_info/output_manifest.tsv').read_text()
            self.assertIn('results/sample.mge_host.filtered.tsv',inventory)
            self.assertNotIn('.snakemake',inventory)
            self.assertTrue((root/'README_outputs.md').is_file())
    def test_declared_offsets(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);(p/'bins').mkdir();(p/'bins/host.fa').write_text('>a\nACG\n>b\nTT\n')
            utils.merge_genomes(p/'bins',str(p/'host_genomes.fasta'),str(p/'host_contig_offsets.tsv'))
            self.assertEqual((p/'host_genomes.fasta').read_text(),'>host\nACGTT\n')
            self.assertIn('host\tb\t3\t5',(p/'host_contig_offsets.tsv').read_text())
