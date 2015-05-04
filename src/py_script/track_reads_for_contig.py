import glob

pread_did_to_rid = open("../1-preads_ovl/preads_ids").read().split("\n")
rid_to_oid = open("../0-rawreads/raw_reads_ids").read().split("\n")

max_overlap_pread = {}
for fn in glob.glob("../1-preads_ovl/max_ovlp.*"):
    with open(fn) as f:
        for row in f:
            row = row.strip().split()
            max_overlap_pread.setdefault(int(row[0]), set())
            max_overlap_pread[int(row[0])].add(int(row[1]))

max_overlap_rread = {}
for fn in glob.glob("../0-rawreads/max_ovlp.*"):
    with open(fn) as f:
        for row in f:
            row = row.strip().split()
            max_overlap_rread.setdefault(int(row[0]), set())
            max_overlap_rread[int(row[0])].add(int(row[1]))

ctg_to_preads = {}
rid_set = set()

for fn in ("p_ctg_tiling_path", "a_ctg_tiling_path"):
    with open(fn) as f:
        for row in f:
            row = row.strip().split()
            ctg = row[0]
            frg0 = int(row[3])
            ctg_to_preads.setdefault( ctg, set() )

            rid = pread_did_to_rid[frg0].split("/")[1]
            rid = int(rid[:-1])
            oid = rid_to_oid[rid]
            ctg_to_preads[ctg].add((0, frg0, rid, oid))
            rid_set.add(frg0)

            for frg in list(max_overlap_pread.get(frg0,[])):
                rid = pread_did_to_rid[frg].split("/")[1]
                rid = int(rid[:-1])
                if rid in rid_set:
                    continue
                oid = rid_to_oid[rid]
                ctg_to_preads[ctg].add((1, frg, rid, oid))
                rid_set.add(rid)


for i in range(2,4):
    for ctg in ctg_to_preads:
        for r in list(ctg_to_preads[ctg]):
            class_, pid, rid0, oid = r

            for rid in list(max_overlap_rread.get(rid0,[])):
                if rid in rid_set:
                    continue
                oid = rid_to_oid[rid]
                ctg_to_preads[ctg].add((i ,rid0, rid, oid))
                rid_set.add(rid)


for ctg in ctg_to_preads:
    for r in list(ctg_to_preads[ctg]):
        print ctg, r[0], "%09d" % r[1], "%09d" % r[2], r[3]
