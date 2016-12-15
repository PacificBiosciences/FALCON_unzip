import pysam

import argparse
import glob
import os
import sys

def select_reads_from_bam(input_bam_fofn_fn, rawread_to_contigs_fn, rawread_ids_fn, sam_dir):
    """Write ctg.sam files into cwd,
    for each 'ctg' read in input BAMs.
    """
    read_partition = {}
    read_to_ctgs = {}

    print "rawread_ids_fn:", repr(rawread_ids_fn)
    print "rawread_to_contigs_fn:", repr(rawread_to_contigs_fn)
    rid_to_oid = open(rawread_ids_fn).read().split('\n')
    with open(rawread_to_contigs_fn) as f:
        for row in f:
            row = row.strip().split()
            if int(row[3]) >= 1: #keep top one hits
                continue
            ctg_id = row[1]
            if ctg_id == 'NA':
                continue
            read_partition.setdefault( ctg_id, set())
            r_id = row[0]
            o_id = rid_to_oid[ int(r_id) ]
            read_partition[ ctg_id ].add(o_id)
            read_to_ctgs.setdefault(o_id, [])
            read_to_ctgs[ o_id ].append( (int(row[4]) ,ctg_id) )
    print "num read_partitions:", len(read_partition)
    print "num read_to_ctgs:", len(read_to_ctgs)

    header = None
    fofn_basedir = os.path.normpath(os.path.dirname(input_bam_fofn_fn))
    def abs_fn(maybe_rel_fn):
        if os.path.isabs(maybe_rel_fn):
            return maybe_rel_fn
        else:
            return os.path.join(fofn_basedir, maybe_rel_fn)
    for row in open(input_bam_fofn_fn):
        fn = abs_fn(row.strip())
        samfile = pysam.AlignmentFile(fn, 'rb', check_sq = False)
        if header == None:
            header = samfile.header
        else:
            header['RG'].extend(samfile.header['RG'])
        samfile.close()

    PG = header.pop('PG') #remove PG line as there might be a bug that generates no readable chrs
    #print PG

    #base_dir = os.getcwd()
    #outfile = pysam.AlignmentFile( os.path.join(base_dir, 'header.sam' ), 'wh', header=header )
    #outfile.close()

    ctgs = read_partition.keys()
    ctgs.sort()
    selected_ctgs = set()
    for ctg in ctgs:
        picked_reads = read_partition[ ctg ]
        print "ctg, len:", ctg, len(picked_reads) # TODO: Is this for debugging?
        if len(picked_reads) > 20:
            selected_ctgs.add(ctg)

    outfile = {}

    for row in open(input_bam_fofn_fn):
        fn = abs_fn(row.strip())
        samfile = pysam.AlignmentFile(fn, 'rb', check_sq = False)
        for r in samfile.fetch( until_eof = True ):
            if r.query_name not in read_to_ctgs:
                #print "Missing:", r.query_name
                continue
            ctg_list = read_to_ctgs[ r.query_name ]
            ctg_list.sort()
            score, ctg = ctg_list[0]
            if ctg not in selected_ctgs:
                #print "Not selected:", ctg
                continue
            if ctg not in outfile:
                samfile_fn = os.path.join(sam_dir, '%s.sam' % ctg)
                print >>sys.stderr, 'samfile_fn:{!r}'.format(samfile_fn)
                outfile[ctg] = pysam.AlignmentFile(samfile_fn, 'wh', header=header)
            outfile[ctg].write(r)
        samfile.close()

    for ctg in outfile:
        outfile[ctg].close()

def parse_args(argv):
    parser = argparse.ArgumentParser(description='Write ctg.sam files, based on BAM subreads.',
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--rawread-to-contigs', type=str, default='./2-asm-falcon/read_maps/rawread_to_contigs', help='rawread_to_contigs file (from where?)')
    parser.add_argument('--rawread-ids', type=str, default='./2-asm-falcon/read_maps/dump_rawread_ids/rawread_ids', help='rawread_ids file (from where?)')
    parser.add_argument('--sam-dir', type=str, default='./4-quiver/reads', help='Output directory for ctg.sam files')
    parser.add_argument('input_bam_fofn', type=str, help='File of BAM filenames. Paths are relative to dir of FOFN, not CWD.')
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)

    select_reads_from_bam(args.input_bam_fofn, args.rawread_to_contigs, args.rawread_ids, args.sam_dir)
