#!/usr/bin/env python3
"""BamM 1.7.3 Python 3 adapter: encode paths passed to its C ABI.

The BamM C filtering algorithm and threshold semantics are unchanged.
No installed third-party file is modified. This module patches the local
module binding only, for the lifetime of this Python process.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile
import bamm.bamFilter as module

original_wrapper = module.CWrapper

def compatible_wrapper():
    wrapper = original_wrapper()
    original_filter = wrapper._filterReads
    def filter_reads(source, destination, *args):
        return original_filter(os.fsencode(source), os.fsencode(destination), *args)
    wrapper._filterReads = filter_reads
    return wrapper

module.CWrapper = compatible_wrapper

def run(source, destination, identity, alignment):
    module.BamFilter(str(source), outFolder=str(destination), minPcId=identity,
                     minPcAln=alignment).filter()

def self_test():
    with tempfile.TemporaryDirectory() as folder:
        p = Path(folder)
        sam = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ref\tLN:10000\n'
        for i, nm in enumerate([0, 1, 2], 1):
            sam += f'read{i}\t0\tref\t{i*100}\t60\t100M\t*\t0\t0\t' + 'ACGT'*25 + '\t' + 'I'*100 + f'\tNM:i:{nm}\n'
        (p/'test.sam').write_text(sam)
        subprocess.run(['samtools','view','-b','-o',str(p/'test.bam'),str(p/'test.sam')],check=True)
        run(p/'test.bam',p,0.99,0.9)
        result = subprocess.check_output(['samtools','view',str(p/'test_filtered.bam')],text=True)
        names = [line.split('\t')[0] for line in result.splitlines()]
        if names != ['read1','read2']:
            raise RuntimeError(f'BamM 99% threshold smoke test failed: {names}')
    print('BamM actual filtering passed: 100% and 99% retained; 98% excluded.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-b')
    parser.add_argument('-o')
    parser.add_argument('--percentage_id',type=float,default=0.99)
    parser.add_argument('--percentage_aln',type=float,default=0.9)
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        self_test()
    else:
        if not args.b or not args.o:
            parser.error('-b and -o required')
        run(args.b,args.o,args.percentage_id,args.percentage_aln)
