#!/usr/bin/env python3
"""
检查 fasta 文件，如果包含多条序列则合并为一条序列输出，否则直接复制
用法: python check_and_merge_fasta.py input.fasta output.fasta
"""

import sys

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

sequences = []
header = None

with open(input_fasta) as f:
    seq_lines = []
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if seq_lines:
                sequences.append("".join(seq_lines))
                seq_lines = []
            header = line[1:].split()[0]
        else:
            seq_lines.append(line)
    if seq_lines:
        sequences.append("".join(seq_lines))

# 如果只有一条序列，直接复制
if len(sequences) == 1:
    with open(output_fasta, "w") as out:
        with open(input_fasta) as f_in:
            out.write(f_in.read())
else:
    # 合并为一条序列
    merged_seq = "".join(sequences)
    base_name = input_fasta.split("/")[-1].split(".")[0]
    with open(output_fasta, "w") as out:
        out.write(f">{base_name}\n{merged_seq}\n")
