import subprocess
import sys
import tempfile
from pathlib import Path
import unittest
import pandas as pd

from metaspici.cli import WORKFLOW
SCRIPT = WORKFLOW / 'scripts/hic_foldchange.py'

class ResultTest(unittest.TestCase):
    def run_case(self, links, cutoff=50, length=1000, read_length=10):
        t = tempfile.TemporaryDirectory()
        self.addCleanup(t.cleanup)
        p = Path(t.name)
        hosts = ['a_host', 'b_host', 'background']
        mges = ['z_mge', 'zz_mge']
        def fasta(ids): return ''.join('>'+i+'\n'+'A'*length+'\n' for i in ids)
        (p/'hosts.fa').write_text(fasta(hosts)); (p/'mges.fa').write_text(fasta(mges))
        (p/'ref.fa').write_text(fasta(hosts+mges))
        (p/'reads.fq').write_text('@x\n'+'A'*read_length+'\n+\n'+'I'*read_length+'\n')
        (p/'links.tsv').write_text(links)
        subprocess.run([sys.executable,str(SCRIPT),str(p/'links.tsv'),str(p/'ref.fa'),
                        str(p/'reads.fq'),str(p/'stats.tsv'),str(p/'final.tsv'),'.14',str(cutoff),'4',
                        str(p/'all.tsv'),str(p/'hosts.fa'),str(p/'mges.fa')],check=True,capture_output=True)
        return pd.read_csv(p/'all.tsv',sep='\t'),pd.read_csv(p/'final.tsv',sep='\t'),pd.read_csv(p/'stats.tsv',sep='\t')

    def test_all_types_and_final_orientation(self):
        links = ('background\tbackground\t1000000\n'
                 'a_host\ta_host\t100\nb_host\tb_host\t100\n'
                 'z_mge\tz_mge\t100\nzz_mge\tzz_mge\t100\n'
                 'a_host\tz_mge\t10\nb_host\tz_mge\t2\n'
                 'a_host\tb_host\t10\nz_mge\tzz_mge\t10\n')
        all_edges, final, _ = self.run_case(links)
        self.assertEqual(len(all_edges),4)
        self.assertTrue((all_edges.contact_value==2).any())
        self.assertEqual(list(zip(final.source,final.target)),[('z_mge','a_host')])
        fe = final.foldchange.iloc[0]
        original = all_edges.loc[(all_edges.source=='a_host') & (all_edges.target=='z_mge'),'foldchange'].iloc[0]
        self.assertEqual(fe, original)
        a2, f2, _ = self.run_case(links,cutoff=fe*1.01)
        self.assertTrue(f2.empty)
        pd.testing.assert_frame_equal(all_edges,a2)

    def test_fixed_depth_and_boundary(self):
        # This test explicitly overrides the default with a 0.14 cutoff.
        # a_host has 14 read ends: 14*10/1000 = 0.14, so it survives.
        # b_host has 13 read ends: 0.13, so it is removed.
        links='a_host\ta_host\t6\na_host\tz_mge\t2\nb_host\tb_host\t6\nb_host\tz_mge\t1\nz_mge\tz_mge\t100\n'
        all_edges, _, stats = self.run_case(links,cutoff=0)
        self.assertEqual(set(all_edges.source),{'a_host'})
        self.assertAlmostEqual(stats.set_index('Genome').loc['a_host','Depth_Total'],.14)
        # At read length 150, the explicitly supplied depth cutoff stays .14, not .02*150.
        a150, _, _ = self.run_case(links,cutoff=0,read_length=150)
        self.assertEqual(len(a150),2)

    def test_empty_network_and_all_filtered(self):
        for links in ('','a_host\tz_mge\t1\n'):
            a,f,_=self.run_case(links)
            self.assertTrue(a.empty and f.empty)
            self.assertEqual(list(f.columns),['source','target','contact_value','foldchange'])

if __name__ == '__main__': unittest.main()
