import os
from falcon_kit.FastaReader import FastaReader

os.system("nucmer -mum p_ctg.fa h_ctg_all.fa -p hp_aln")
os.system("show-coords -T -H -l -c hp_aln.delta > hp_aln.coor")

filter_out = set()
with open("hp_aln.coor") as f:
    for row in f:
        row = row.strip().split()
        q_cov = float(row[10])
        idt = float(row[6])
        if q_cov > 99 and idt > 99:
            filter_out.add(row[-1])

with open("h_ctg.fa", "w") as f:
    h_tig_all = FastaReader("h_ctg_all.fa")
    for r in h_tig_all:
        if r.name in filter_out:
            continue
        print >>f, ">"+r.name
        print >>f, r.sequence
