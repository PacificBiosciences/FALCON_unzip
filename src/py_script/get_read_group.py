
rid_map = {}
with open("q_id_map") as f:
    for l in f:
        l = l.strip().split()
        rid_map[int(l[0])] = l[1]


read_to_variants = {}
variant_to_reads = {}
with open("vmap2") as f: 
    for l in f:
        l = l.strip().split()
        variant = "_".join(l[:3])
        read_id = int(l[3])
        read_to_variants.setdefault(read_id, set())
        read_to_variants[read_id].add(variant)
        variant_to_reads.setdefault(variant, set())
        variant_to_reads[variant].add(read_id)


variant_to_phase = {}
with open("phased_variants") as f:
    for l in f:
        l = l.strip().split()
        """V 1 6854 6854_A_A 6854_A_G 6854 22781"""
        if l[0] != "V":
            continue
        pb_id = int(l[1])
        variant_to_phase[ l[3] ] = (pb_id, 0)
        variant_to_phase[ l[4] ] = (pb_id, 1)

for r in read_to_variants:
    vl = {}
    pl = set()
    for v in list( read_to_variants[r] ):
        if v in variant_to_phase:
            p = variant_to_phase[v]
            vl[ p ] = vl.get(p, 0) + 1
            pl.add(p[0])
    pl = list(pl)
    pl.sort()
    for p in pl:
        if vl.get( (p,0), 0) - vl.get( (p,1), 0) > 1:
            print r, p, 0, vl.get( (p,0), 0), vl.get( (p,1), 0), rid_map[r]
        elif vl.get( (p,1), 0) - vl.get( (p,0), 0) > 1:
            print r, p, 1, vl.get( (p,0), 0), vl.get( (p,1), 0), rid_map[r]
        


