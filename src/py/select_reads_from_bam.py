import pysam
import sys
import glob
import os

def main(args):
    read_partition = {}
    read_to_ctgs = {}

    fn = "./3-unzip/read_maps/rawread_to_contigs"
    rid_to_oid = open("2-asm-falcon/read_maps/raw_read_ids").read().split("\n")
    with open(fn) as f:
        for row in f:
            row = row.strip().split()
            if int(row[3]) >= 1: #keep top one hits
                continue
            ctg_id = row[1]
            if ctg_id == "NA":
                continue
            read_partition.setdefault( ctg_id, set() )
            r_id = row[0]
            o_id = rid_to_oid[ int(r_id) ]
            read_partition[ ctg_id ].add( o_id )
            read_to_ctgs.setdefault( o_id, [] )
            read_to_ctgs[ o_id ].append( (int(row[4]) ,ctg_id) )

    header = None
    for row in open(sys.argv[1]):
        fn = row.strip()
        samfile = pysam.AlignmentFile(fn, "rb", check_sq = False )
        if header == None: 
            header = samfile.header
        else:
            header["RG"].extend( samfile.header["RG"] )
        samfile.close()

    PG = header.pop("PG") #remove PG line as there might be a bug that generates no readable chrs
    #print PG 

    base_dir = os.getcwd()
    #outfile = pysam.AlignmentFile( os.path.join(base_dir, "header.sam" ), "wh", header=header )
    #outfile.close()

    ctgs = read_partition.keys()
    ctgs.sort()
    selected_ctgs = set()
    for ctg in ctgs:
        picked_reads = read_partition[ ctg ]
        print ctg, len(picked_reads)
        if len(picked_reads) > 20:
            selected_ctgs.add(ctg)

    outfile = {}

    for row in open(sys.argv[1]):
        fn = row.strip()
        samfile = pysam.AlignmentFile(fn, "rb", check_sq = False )
        for r in samfile.fetch( until_eof = True ):
            if r.query_name not in read_to_ctgs:
                continue
            ctg_list = read_to_ctgs[ r.query_name ]
            ctg_list.sort()
            score, ctg = ctg_list[0]
            if ctg not in selected_ctgs:
                continue
            if ctg not in outfile:
                outfile[ctg] = pysam.AlignmentFile( os.path.join(base_dir, "4-quiver/reads/", "%s.sam" % ctg), "wh", header=header )
            outfile[ctg].write(r)
        samfile.close()

    for ctg in outfile:
        outfile[ctg].close()
