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

p_ctg_to_phase = {}
with open("p_ctg_path") as f:
    for row in f:
        row = row.strip().split()
        b_id, ph_id = (int(row[-2]), int(row[-1]) )
        p_ctg_to_phase.setdefault( row[0], {} )
        p_ctg_to_phase[row[0]].setdefault( ( b_id, ph_id ), 0)
        p_ctg_to_phase[row[0]][ ( b_id, ph_id ) ] += 1


h_ctg_to_phase = {}
with open("h_ctg_path") as f:
    for row in f:
        row = row.strip().split()
        b_id, ph_id = (int(row[-2]), int(row[-1]) )
        h_ctg_to_phase.setdefault( row[0], {} )
        h_ctg_to_phase[row[0]].setdefault( ( b_id, ph_id ), 0)
        h_ctg_to_phase[row[0]][ ( b_id, ph_id ) ] += 1

h_ids = open("h_ctg_ids","w")
with open("h_ctg.fa", "w") as f:
    h_tig_all = FastaReader("h_ctg_all.fa")
    for r in h_tig_all:
        p_ctg_phase = p_ctg_to_phase.get(r.name.split("_")[0], {})

        if r.name in filter_out:
            edge_count = sum(h_ctg_to_phase[ r.name ].values())
            same_phase_to_p_ctg_count = 0
            for  b_id, ph_id in h_ctg_to_phase[ r.name ]:
                if b_id != -1:
                    if (b_id, ph_id) in p_ctg_phase:
                        same_phase_to_p_ctg_count += h_ctg_to_phase[ r.name ][ (b_id, ph_id) ]
            unphased_edge_count = h_ctg_to_phase[ r.name ] .get( (-1, 0), 0 )

            print r.name, edge_count, unphased_edge_count, same_phase_to_p_ctg_count
            if edge_count - unphased_edge_count - same_phase_to_p_ctg_count < 5: # there are many non-p_ctg phase segment, do not filter out
                continue

        print >>f, ">"+r.name
        print >>f, r.sequence
        print >> h_ids, r.name
h_ids.close()
