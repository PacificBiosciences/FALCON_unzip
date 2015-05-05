import os

if __name__ == "__main__":
    import argparse
    import re
    parser = argparse.ArgumentParser(description='mapping internal daligner read id to phase block and phase')
    # we can run this in parallel mode in the furture
    #parser.add_argument('--n_core', type=int, default=4,
    #                    help='number of processes used for generating consensus')
    parser.add_argument('--phased_reads', type=str, help='path to read vs. phase map', required=True)
    parser.add_argument('--ctg_map', type=str, help='path to the read-contig map files', required=True) 
    parser.add_argument('--ctg_id', type=str, help='contig identifier in the bam file', required=True)
    parser.add_argument('--base_dir', type=str, default="./", help='the output base_dir, default to current working directory')

    args = parser.parse_args()
    phased_reads = args.phased_reads
    ctg_map = args.ctg_map
    ctg_id = args.ctg_id
    base_dir = args.base_dir

    rid_to_phase = {}        
    with open(phased_reads) as f:
        for row in f:
            row = row.strip().split()
            rid_to_phase[row[5]] = (int(row[1]), int(row[2]))

    arid_to_phase = {}        
    with open(ctg_map) as f:
        for row in f:
            row = row.strip().split()
            ctg_id = row[0]
            if not ctg_id.startswith(ctg_id):
                continue
            if int(row[1]) > 1: #not preads
                continue
            phase = rid_to_phase.get( row[4], (-1, 0) )
            arid_to_phase["%09d" % int(row[2])] = phase
            
    with open(os.path.join(base_dir, "rid_to_phase"),"w") as f:
        for arid, phase in arid_to_phase.items():
            print >>f, arid, phase[0], phase[1]
