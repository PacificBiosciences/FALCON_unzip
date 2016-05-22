from falcon_kit.FastaReader import FastaReader
import os
import glob
import sys


rawread_dir = os.path.abspath( "./0-rawreads" )
pread_dir = os.path.abspath( "./1-preads_ovl" )
asm_dir = os.path.abspath( "./2-asm-falcon" )

def main():
    import argparse
    import re
    parser = argparse.ArgumentParser(description='align reads for one contig')
    # we can run this in parallel mode in the furture
    #parser.add_argument('--n_core', type=int, default=4,
    #                    help='number of processes used for generating consensus')
    parser.add_argument('--fofn', type=str, default="./input.fofn", help='path to the file of the list of raw read fasta files')
    parser.add_argument('--ctg_fa', type=str, default="./2-asm-falcon/p_ctg.fa", help='path to the contig fasta file')
    parser.add_argument('--ctg_id', type=str, default="all", help='contig identifier in the contig fasta file')
    parser.add_argument('--read_map_dir', type=str, default="./2-asm-falcon/read_maps", help='path to the read-contig map directory') 
    parser.add_argument('--base_dir', default="./3-unzip/reads", type=str, help='the output base_dir, default to current working directory')

    rawread_id_file = os.path.join( rawread_dir, "raw_read_ids" )
    pread_id_file = os.path.join( pread_dir, "pread_ids" ) 

    rid_to_oid = open(rawread_id_file).read().split("\n")  #daligner raw read id to the original ids
    pid_to_fid = open(pread_id_file).read().split("\n")  #daligner pread id to the fake ids

    def pid_to_oid(pid):
        fid = pid_to_fid[int(pid)]
        rid = int(fid.split("/")[1])/10
        return rid_to_oid[int(rid)]
    
    
    
    args = parser.parse_args()
    read_fofn = args.fofn
    ctg_fa = args.ctg_fa
    ctg_id = args.ctg_id
    read_map_dir = args.read_map_dir
    base_dir = args.base_dir

    ref_fasta = FastaReader(ctg_fa)
    all_ctg_ids = set()
    for s in ref_fasta:
        s_id = s.name.split()[0]
        if ctg_id != "all" and s_id != ctg_id:
            continue

        if len(s.sequence) < 20000:
            continue
        if ctg_id != "all":
            ref_out = open( os.path.join( base_dir, "%s_ref.fa" % ctg_id), "w" )
        else:
            ref_out = open( os.path.join( base_dir, "%s_ref.fa" % s_id), "w" )

        print >>ref_out, ">%s" % s_id
        print >>ref_out, s.sequence
        all_ctg_ids.add(s_id)
        ref_out.close()

    
    read_set = {}
    ctg_id_hits = {} 

    map_fn = os.path.join(read_map_dir,"rawread_to_contigs")
    with open(map_fn, "r") as f:
        for row in f:
            row = row.strip().split()
            hit_ctg = row[1]
            hit_ctg = hit_ctg.split("-")[0]
            if int(row[3]) == 0:
                o_id = rid_to_oid[int(row[0])]
                read_set[o_id] = hit_ctg
                ctg_id_hits[hit_ctg] = ctg_id_hits.get(hit_ctg, 0) + 1

    map_fn = os.path.join(read_map_dir,"pread_to_contigs")
    with open(map_fn, "r") as f:
        for row in f:
            row = row.strip().split()
            hit_ctg = row[1]
            hit_ctg = hit_ctg.split("-")[0]
            if hit_ctg not in read_set and int(row[3]) == 0:
                o_id = pid_to_oid(row[0])
                read_set[o_id] = hit_ctg
                ctg_id_hits[hit_ctg] = ctg_id_hits.get(hit_ctg, 0) + 1
                

    with open(os.path.join( base_dir, "ctg_list"),"w") as f:
        for ctg_id in sorted(list(all_ctg_ids)):
            if ctg_id_hits.get(ctg_id, 0) < 5:
                continue
            if ctg_id[-1] not in ["F", "R"]: #ignore small circle contigs, they need different approach
                continue
            print >>f, ctg_id

    read_out_files = {}
    with open(read_fofn, "r") as f:
        for r_fn in f:
            r_fn = r_fn.strip()
            read_fa_file = FastaReader(r_fn)
            for r in read_fa_file:
                rid = r.name.split()[0]
                if rid not in read_set:
                    ctg_id = "unassigned"
                else:     
                    ctg_id = read_set[rid]

                if ctg_id == "NA" or ctg_id not in all_ctg_ids:
                    ctg_id = "unassigned"

                if ctg_id not in read_out_files:
                    read_out = open( os.path.join( base_dir, "%s_reads.fa" % ctg_id), "w" )
                    read_out_files[ctg_id] = read_out
                else:
                    read_out = read_out_files[ctg_id]

                print >>read_out, ">"+rid
                print >>read_out, r.sequence
    
    for read_out in read_out_files.values():
        read_out.close()
