import sys
import gzip

fq = sys.argv[1]

def open_file(f):
    return gzip.open(f, "rt") if f.endswith(".gz") else open(f)

count = 0
total = 0
max_reads = 100000

with open_file(fq) as f:
    for i, line in enumerate(f):
        if i % 4 == 1:
            total += len(line.strip())
            count += 1
            if count >= max_reads:
                break

mean_len = total / count
cutoff = int(mean_len * 0.8)

print(cutoff)
