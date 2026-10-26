import sys
import pandas as pd
import numpy as np
from collections import defaultdict
import itertools
from Bio import SeqIO

# ===== 输入参数 =====
hic_file = sys.argv[1]
fasta_file = sys.argv[2]
fastq_file = sys.argv[3]
stats_out = sys.argv[4]
result_out = sys.argv[5]
all_out = sys.argv[9]
host_ids = {r.id for r in SeqIO.parse(sys.argv[10], 'fasta')}
mge_ids = {r.id for r in SeqIO.parse(sys.argv[11], 'fasta')}
if host_ids & mge_ids:
    raise ValueError('Host and MGE reference IDs overlap')

def write_empty_results():
    for path in (all_out, result_out):
        pd.DataFrame(columns=['source','target','contact_value','foldchange']).to_csv(path, sep='\t', index=False)

# ===== Step 1. 读取 Hi-C =====
df = pd.read_csv(hic_file, sep='\t',
                 names=['source', 'target', 'contact_value'])

# ===== Step 2. 估算 read length =====
def estimate_read_length(fastq_file, max_reads=10000):
    lengths = []
    with open(fastq_file) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                lengths.append(len(line.strip()))
                if len(lengths) >= max_reads:
                    break
    return int(np.mean(lengths)) if lengths else 0

read_length = estimate_read_length(fastq_file)
depth_threshold = float(sys.argv[6])

print("Read length:", read_length)
print("Depth threshold:", depth_threshold)

# ===== Step 3. genome-level reads =====
total_reads = {}
self_reads = {}

for _, row in df.iterrows():
    g1 = row['source']
    g2 = row['target']
    pairs = row['contact_value']

    if g1 == g2:
        reads = pairs * 2
        total_reads[g1] = total_reads.get(g1, 0) + reads
        self_reads[g1] = self_reads.get(g1, 0) + reads
    else:
        total_reads[g1] = total_reads.get(g1, 0) + pairs
        total_reads[g2] = total_reads.get(g2, 0) + pairs

genomes = list(total_reads.keys())

stats_df = pd.DataFrame({
    'Genome': genomes,
    'TotalReads': [total_reads[g] for g in genomes],
    'SelfReads': [self_reads.get(g, 0) for g in genomes],
})

stats_df['NonSelfReads'] = stats_df['TotalReads'] - stats_df['SelfReads']

if df.empty:
    pd.DataFrame(columns=['Genome','TotalReads','SelfReads','NonSelfReads','GenomeLength','Depth_Total']).to_csv(stats_out, sep='\t', index=False)
    write_empty_results()
    sys.exit(0)

# ===== Step 4. genome length =====
genome_lengths = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    genome_lengths[record.id] = len(record.seq)

length_df = pd.DataFrame({
    'Genome': list(genome_lengths.keys()),
    'GenomeLength': list(genome_lengths.values())
})

# ===== Step 5. depth =====
merged_stats = stats_df.merge(length_df, on='Genome', how='left')

merged_stats['Depth_Total'] = (
    merged_stats['TotalReads'] * read_length / merged_stats['GenomeLength']
)

merged_stats.to_csv(stats_out, sep='\t', index=False)

# ===== Step 6. 过滤低丰度 =====
low_depth = merged_stats.loc[
    merged_stats['Depth_Total'] < depth_threshold,
    'Genome'
].tolist()

df = df[
    ~(df['source'].isin(low_depth) |
      df['target'].isin(low_depth))
].copy()

# Empty input/all genomes filtered: valid empty result with stable headers.
if df.empty:
    write_empty_results()
    sys.exit(0)

# ===== Step 7. masked ratio =====
masked_ratio = {}
for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq).upper()
    masked_ratio[record.id] = seq.count('N') / len(seq)

# ===== Step 8. abundance =====
abundance = defaultdict(float)

for _, row in df.iterrows():
    s = row['source']
    t = row['target']
    pairs = row['contact_value']

    if s == t:
        rs = pairs * 2
        rt = 0
    else:
        rs = pairs
        rt = pairs

    fs = masked_ratio.get(s, 0)
    ft = masked_ratio.get(t, 0)

    rs = rs / (1 - fs) if fs < 1 else np.nan
    if rt > 0:
        rt = rt / (1 - ft) if ft < 1 else np.nan

    abundance[s] += rs
    if rt > 0:
        abundance[t] += rt

abundance_df = pd.DataFrame.from_dict(
    abundance, orient='index', columns=['abundance']
)

contigs = abundance_df.index.tolist()

# ===== Step 9. expected =====
pairs_list = []
for g1, g2 in itertools.combinations_with_replacement(contigs, 2):
    a1 = abundance_df.loc[g1, 'abundance']
    a2 = abundance_df.loc[g2, 'abundance']
    pairs_list.append((tuple(sorted([g1, g2])), a1 * a2))

prob_df = pd.DataFrame(pairs_list, columns=['pair', 'raw_prob'])
prob_df['normalized_prob'] = prob_df['raw_prob'] / prob_df['raw_prob'].sum()

# ===== Step 10. observed =====
df['pair'] = df.apply(lambda r: tuple(sorted([r['source'], r['target']])), axis=1)

grouped = df.groupby('pair', as_index=False)['contact_value'].sum()
grouped[['source', 'target']] = pd.DataFrame(grouped['pair'].tolist(), index=grouped.index)

grouped['reads'] = grouped.apply(
    lambda r: r['contact_value'] * 2 if r['source'] == r['target'] else r['contact_value'],
    axis=1
)

merged = grouped.merge(prob_df, on='pair')

total_contacts = merged[merged['source'] != merged['target']]['reads'].sum()

merged['expected'] = merged['normalized_prob'] * total_contacts
merged['foldchange'] = merged['reads'] / merged['expected']

# ===== Step 11. Two association tables =====
# Preserve FE calculation over the complete depth-filtered network above.
# The unfiltered table retains all non-self edges, including 1-3 contacts.
all_edges = merged.loc[merged['source'] != merged['target'],
                       ['source', 'target', 'reads', 'foldchange']].copy()
all_edges = all_edges.rename(columns={'reads': 'contact_value'})
all_edges = all_edges.sort_values('foldchange', ascending=False)
all_edges.to_csv(all_out, sep='\t', index=False)

is_mge_host = (
    (all_edges['source'].isin(mge_ids) & all_edges['target'].isin(host_ids)) |
    (all_edges['source'].isin(host_ids) & all_edges['target'].isin(mge_ids))
)
final = all_edges.loc[
    is_mge_host &
    (all_edges['contact_value'] >= int(sys.argv[8])) &
    (all_edges['foldchange'] > float(sys.argv[7]))
].copy()
# Consistent orientation: source=MGE and target=host, using actual reference IDs.
swap = final['source'].isin(host_ids)
final.loc[swap, ['source', 'target']] = final.loc[swap, ['target', 'source']].to_numpy()
final.to_csv(result_out, sep='\t', index=False)
