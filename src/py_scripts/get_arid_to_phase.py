rid_to_phase = {}        
with open("../4-contig-phasing/phased_reads") as f:
    for row in f:
        row = row.strip().split()
        rid_to_phase[row[5]] = (int(row[1]), int(row[2]))

arid_to_phase = {}        
with open("contig_to_read_map") as f:
    for row in f:
        row = row.strip().split()
        ctg_id = row[0]
        if not ctg_id.startswith("000000F"):
            continue
        if int(row[1]) > 1: #not preads
            continue
        phase = rid_to_phase.get( row[4], (-1, 0) )
        arid_to_phase["%09d" % int(row[2])] = phase
        
with open("arid_to_phase","w") as f:
    for arid, phase in arid_to_phase.items():
        print >>f, arid, phase[0], phase[1]
