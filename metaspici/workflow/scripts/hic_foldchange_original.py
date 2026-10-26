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
depth_threshold = 0.02 * read_length

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

# ===== Step 11. 输出 =====
# ===== Step 11. 输出（筛选 foldchange > 30）=====
final = merged[
    (merged['source'] != merged['target']) &
    (merged['reads'] > 3)
].copy()

# 👉 过滤 foldchange > 30
final = final[final['foldchange'] > 30]

# 👉 只保留需要的列
final = final[['source', 'target', 'reads', 'foldchange']]

# 👉 排序
final = final.sort_values('foldchange', ascending=False)

# 👉 重命名 reads → contact_value（更符合你输入格式）
final = final.rename(columns={'reads': 'contact_value'})

# 👉 输出
final.to_csv(result_out, sep='\t', index=False)
