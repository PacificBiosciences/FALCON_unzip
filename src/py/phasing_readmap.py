import os
import sys
import re
import argparse


def get_phasing_readmap(args):
    
    phased_reads = args.phased_reads
    read_map_dir = args.read_map_dir
    the_ctg_id = args.ctg_id
    base_dir = args.base_dir
    
    rawread_id_file = os.path.join( read_map_dir, "raw_read_ids" )
    pread_id_file = os.path.join( read_map_dir, "pread_ids" ) 
    rid_to_oid = open(rawread_id_file).read().split("\n")  #daligner raw read id to the original ids
    pid_to_fid = open(pread_id_file).read().split("\n")  #daligner pread id to the fake ids
    def pid_to_oid(pid):
        fid = pid_to_fid[int(pid)]
        rid = int(fid.split("/")[1])/10
        return rid_to_oid[int(rid)]

    rid_to_oid = open(rawread_id_file).read().split("\n")  #daligner raw read id to the original ids
    pid_to_fid = open(pread_id_file).read().split("\n")  #daligner pread id to the fake ids


    rid_to_phase = {}        
    with open(phased_reads) as f:
        for row in f:
            row = row.strip().split()
            rid_to_phase[row[6]] = (int(row[2]), int(row[3]))

    arid_to_phase = {}        
    map_fn = os.path.join(read_map_dir, "pread_to_contigs")
    with open(map_fn) as f:
        for row in f:
            row = row.strip().split()
            ctg_id = row[1]
            if not ctg_id.startswith(the_ctg_id):
                continue
            if int(row[3]) != 0: #not the best hit
                continue
            o_id = pid_to_oid(row[0])
            phase = rid_to_phase.get( o_id, (-1, 0) )
            arid_to_phase["%09d" % int(row[0])] = phase
            
    with open(os.path.join(base_dir, "rid_to_phase.%s" % the_ctg_id),"w") as f:
        for arid, phase in arid_to_phase.items():
            print >>f, arid, the_ctg_id, phase[0], phase[1]



def parse_args(argv):
    parser = argparse.ArgumentParser(description='mapping internal daligner read id to phase block and phase',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # we can run this in parallel mode in the furture
    #parser.add_argument('--n_core', type=int, default=4,
    #                    help='number of processes used for generating consensus')
    parser.add_argument('--phased_reads', type=str, help='path to read vs. phase map', required=True)
    parser.add_argument('--read_map_dir', type=str, help='path to the read map directory', required=True) 
    parser.add_argument('--ctg_id', type=str, help='contig identifier in the bam file', required=True)
    parser.add_argument('--base_dir', type=str, default="./", help='the output base_dir, default to current working directory')

    args = parser.parse_args(argv[1:])
    return args
    

def main(argv=sys.argv):
    args = parse_args(argv)
    get_phasing_readmap(args)
