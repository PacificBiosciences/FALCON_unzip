import os
import glob

if __name__ == "__main__":
    import argparse
    import re
    parser = argparse.ArgumentParser(description='mapping internal daligner read id to phase block and phase')
    # we can run this in parallel mode in the furture
    #parser.add_argument('--n_core', type=int, default=4,
    #                    help='number of processes used for generating consensus')
    parser.add_argument('--phased_reads', type=str, help='path to read vs. phase map', required=True)
    parser.add_argument('--read_map_dir', type=str, help='path to the read map directory', required=True) 
    parser.add_argument('--ctg_id', type=str, help='contig identifier in the bam file', required=True)
    parser.add_argument('--base_dir', type=str, default="./", help='the output base_dir, default to current working directory')

    args = parser.parse_args()
    phased_reads = args.phased_reads
    read_map_dir = args.read_map_dir
    the_ctg_id = args.ctg_id
    base_dir = args.base_dir

    rid_to_phase = {}        
    with open(phased_reads) as f:
        for row in f:
            row = row.strip().split()
            rid_to_phase[row[6]] = (int(row[2]), int(row[3]))

    arid_to_phase = {}        
    for map_fn in glob.glob(os.path.join(read_map_dir,"pread_to_contigs.*")):
        with open(map_fn) as f:
            for row in f:
                row = row.strip().split()
                ctg_id = row[2]
                if not ctg_id.startswith(the_ctg_id):
                    continue
                if int(row[4]) != 0: #not the best hit
                    continue
                phase = rid_to_phase.get( row[1], (-1, 0) )
                arid_to_phase["%09d" % int(row[0])] = phase
            
    with open(os.path.join(base_dir, "rid_to_phase"),"w") as f:
        for arid, phase in arid_to_phase.items():
            print >>f, arid, phase[0], phase[1]
