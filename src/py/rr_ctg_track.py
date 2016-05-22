from __future__ import print_function
from falcon_kit.multiproc import Pool
import falcon_kit.util.io as io
import argparse
import sys
import glob
import os
from heapq import heappush, heappop, heappushpop

Reader = io.CapturedProcessReaderContext

rawread_dir = os.path.abspath( "./0-rawreads" )
pread_dir = os.path.abspath( "./1-preads_ovl" )
asm_dir = os.path.abspath( "./2-asm-falcon" )

def get_rid_to_ctg(fn):
    rid_to_ctg = {}
    with open(fn) as f:
        for row in f:
            row = row.strip().split()
            pid, rid, oid, ctg = row
            rid_to_ctg.setdefault( rid, set()  )
            rid_to_ctg[ rid ].add(ctg)
    return rid_to_ctg

def run_tr_stage1(db_fn, fn, min_len, bestn, rid_to_ctg):
    cmd = "LA4Falcon -m %s %s" % (db_fn, fn)
    reader = Reader(cmd)
    with reader:
        return fn, tr_stage1(reader.readlines, min_len, bestn, rid_to_ctg)

def tr_stage1(readlines, min_len, bestn, rid_to_ctg):
    rtn = {}
    for l in readlines():
        l = l.strip().split()
        q_id, t_id = l[:2]
        overlap_len = -int(l[2])
        idt = float(l[3])
        q_s, q_e, q_l = int(l[5]), int(l[6]), int(l[7])
        t_s, t_e, t_l = int(l[9]), int(l[10]), int(l[11])
        if t_l < min_len:
            continue
        if q_id not in rid_to_ctg:
            continue
        rtn.setdefault(t_id, [])
        if len(rtn[t_id]) < bestn:
            heappush(rtn[t_id], (overlap_len, q_id) )
        else:
            heappushpop(rtn[t_id], (overlap_len, q_id) )

    return rtn

def run_track_reads(exe_pool, file_list, min_len, bestn, db_fn):
    io.LOG('preparing tr_stage1')
    io.logstats()
    rid_to_ctg = get_rid_to_ctg( os.path.join(asm_dir, "read_maps/read_to_contig_map" ) )
    inputs = []
    for fn in file_list:
        if len(fn) != 0:
            inputs.append( (run_tr_stage1, db_fn, fn, min_len, bestn, rid_to_ctg) )

    bread_to_areads = {}
    for fn, res in exe_pool.imap(io.run_func, inputs):
        for k in res:
            bread_to_areads.setdefault(k, [])
            for item in res[k]:
                if len(bread_to_areads[k]) < bestn:
                    heappush( bread_to_areads[k], item )
                else:
                    heappushpop( bread_to_areads[k], item )

    #rid_to_oid = open(os.path.join( rawread_dir, "raw_read_ids" ) ).read().split("\n")

    with open( os.path.join(asm_dir, "read_maps/rawread_to_contigs"), "w") as out_f:
        for bread in bread_to_areads:

            ctg_score = {}
            for s, rid in bread_to_areads[bread]:
                if rid not in rid_to_ctg: continue
                
                ctgs = rid_to_ctg[rid]
                for ctg in ctgs:
                    ctg_score.setdefault(ctg, [0,0])
                    ctg_score[ctg][0] += -s
                    ctg_score[ctg][1] += 1

            #oid = rid_to_oid[int(bread)]
            ctg_score = ctg_score.items()
            ctg_score.sort( key = lambda k: k[1][0] )
            rank = 0

            for ctg, score_count in ctg_score:
                if bread in rid_to_ctg and ctg in rid_to_ctg[bread]:
                    in_ctg = 1
                else:
                    in_ctg = 0
                score, count = score_count
                #print(bread, oid, ctg, count, rank, score, in_ctg, file=out_f) 
                print(bread, ctg, count, rank, score, in_ctg, file=out_f) 
                rank += 1

        

def try_run_track_reads(n_core, min_len, bestn):
    io.LOG('starting track_reads')
    file_list = glob.glob( os.path.join( rawread_dir, "m*/raw_reads.*.las") )
    io.LOG('file list: %r' % file_list)
    db_fn = os.path.join( rawread_dir, "raw_reads.db" )
    n_core = min(n_core, len(file_list))
    exe_pool = Pool(n_core)
    try:
        run_track_reads(exe_pool, file_list, min_len, bestn, db_fn)
        io.LOG('finished track_reads')
    except KeyboardInterrupt:
        io.LOG('terminating track_reads workers...')
        exe_pool.terminate()

def track_reads(n_core, min_len, bestn, debug, silent, stream):
    if debug:
        n_core = 0
        silent = False
    if silent:
        io.LOG = io.write_nothing
    if stream:
        global Reader
        Reader = io.StreamedProcessReaderContext
    try_run_track_reads(n_core, min_len, bestn)

def parse_args(argv):
    parser = argparse.ArgumentParser(description='a simple multi-processes LAS ovelap data filter',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--n_core', type=int, default=4,
                        help='number of processes used for generating consensus; '
                        '0 for main process only')
    #parser.add_argument('--fofn', type=str, help='file contains the path of all LAS file to be processed in parallel')
    #parser.add_argument('--db', type=str, dest='db_fn', help='read db file path')
    parser.add_argument('--min_len', type=int, default=2500, help="min length of the reads")
    parser.add_argument('--stream', action='store_true', help='stream from LA4Falcon, instead of slurping all at once; can save memory for large data')
    parser.add_argument('--debug', '-g', action='store_true', help="single-threaded, plus other aids to debugging")
    parser.add_argument('--silent', action='store_true', help="suppress cmd reporting on stderr")
    parser.add_argument('--bestn', type=int, default=40, help="keep best n hits")
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    track_reads(**vars(args))

