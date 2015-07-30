from falcon_kit.FastaReader import FastaReader
import os
import glob
import sys

if __name__ == "__main__":
    import argparse
    import re
    parser = argparse.ArgumentParser(description='align reads for one contig')
    # we can run this in parallel mode in the furture
    #parser.add_argument('--n_core', type=int, default=4,
    #                    help='number of processes used for generating consensus')
    parser.add_argument('--fofn', type=str, help='path to the file of the list of raw read fasta files', required=True)
    parser.add_argument('--ctg_fa', type=str, help='path to the contig fasta file', required=True)
    parser.add_argument('--ctg_id', type=str, help='contig identifier in the contig fasta file', required=True)
    parser.add_argument('--read_map_dir', type=str, help='path to the read-contig map directory', required=True) 
    parser.add_argument('--base_dir', type=str, default="./", help='the output base_dir, default to current working directory')

    args = parser.parse_args()
    read_fofn = args.fofn
    ctg_fa = args.ctg_fa
    ctg_id = args.ctg_id
    read_map_dir = args.read_map_dir
    base_dir = args.base_dir

    read_set = set()
    for map_fn in glob.glob(os.path.join(read_map_dir,"rawread_to_contigs.*")):
        with open(map_fn, "r") as f:
            for row in f:
                row = row.strip().split()
                if row[2].startswith(ctg_id) and int(row[4]) == 0:
                    read_set.add(row[1])
                
    ref_fasta = FastaReader(ctg_fa)
    ref_out = open( os.path.join( base_dir, "%s_ref.fa" % ctg_id), "w" )
    for s in ref_fasta:
        s_id = s.name.split()[0]
        if s_id != ctg_id:
            continue
        print >>ref_out, ">%s" % ctg_id
        print >>ref_out, s.sequence
    ref_out.close()


    read_out = open( os.path.join( base_dir, "%s_reads.fa" % ctg_id), "w" )
    with open(read_fofn, "r") as f:
        for r_fn in f:
            r_fn = r_fn.strip()
            read_fa_file = FastaReader(r_fn)
            for r in read_fa_file:
                rid = r.name.split()[0]
                if rid not in read_set:
                    continue
                print >>read_out, ">"+rid
                print >>read_out, r.sequence

    read_out.close()
