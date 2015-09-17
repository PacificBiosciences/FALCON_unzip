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


h_ctg_to_phase = {}
with open("h_ctg_path") as f:
    for row in f:
        row = row.strip().split()
        b_id, ph_is = (int(row[-2]), int(row[-1]) )
        h_ctg_to_phase.setdefault( row[0], {} )
        h_ctg_to_phase[row[0]].setdefault( ( b_id, ph_is ), 0)
        h_ctg_to_phase[row[0]][ ( b_id, ph_is ) ] += 1

h_ids = open("h_ctg_ids","w")
with open("h_ctg.fa", "w") as f:
    h_tig_all = FastaReader("h_ctg_all.fa")
    for r in h_tig_all:
        if r.name in filter_out:
            edge_count = sum(h_ctg_to_phase[ r.name ].values())
            unphased_edge_count = h_ctg_to_phase[ r.name ] .get( (-1, 0), 0 )

            print r.name, edge_count, unphased_edge_count
            if edge_count - unphased_edge_count < 5:
                continue

        print >>f, ">"+r.name
        print >>f, r.sequence
        print >> h_ids, r.name
h_ids.close()
