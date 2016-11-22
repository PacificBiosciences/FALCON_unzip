#from pypeflow.common import *
#from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
#from pypeflow.task import PypeTask, PypeThreadTaskBase, PypeTaskBase
#from pypeflow.controller import PypeWorkflow, PypeMPWorkflow, PypeThreadWorkflow
from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase,
        makePypeLocalFile, fn, PypeTask)
PypeThreadTaskBase = MyFakePypeThreadTaskBase
from falcon_kit.FastaReader import FastaReader
from falcon_kit.fc_asm_graph import AsmGraph
import glob
import sys
import subprocess as sp
import shlex
import os

def make_dirs(d):
    if not os.path.isdir(d):
        os.makedirs(d)

rawread_dir = os.path.abspath( "./0-rawreads" )
pread_dir = os.path.abspath( "./1-preads_ovl" )
asm_dir = os.path.abspath( os.path.join("./3-unzip/") )

read_map_dir = os.path.abspath(os.path.join(asm_dir, "read_maps"))
make_dirs(read_map_dir)

wf = PypeProcWatcherWorkflow(
        max_jobs=12,
)

rawread_db = makePypeLocalFile( os.path.join( rawread_dir, "raw_reads.db" ) )
rawread_id_file = makePypeLocalFile( os.path.join( rawread_dir, "raw_read_ids" ) )

@PypeTask( inputs = {"rawread_db": rawread_db}, 
           outputs =  {"rawread_id_file": rawread_id_file},
           TaskType = PypeThreadTaskBase,
           URL = "task://localhost/dump_rawread_ids" )
def dump_rawread_ids(self):
    rawread_db = fn( self.rawread_db )
    rawread_id_file = fn( self.rawread_id_file )
    os.system("DBshow -n %s | tr -d '>' | awk '{print $1}' > %s" % (rawread_db, rawread_id_file) )

wf.addTask( dump_rawread_ids )

pread_db = makePypeLocalFile( os.path.join( pread_dir, "preads.db" ) )
pread_id_file = makePypeLocalFile( os.path.join( pread_dir, "pread_ids" ) )

@PypeTask( inputs = {"pread_db": pread_db}, 
           outputs =  {"pread_id_file": pread_id_file},
           TaskType = PypeThreadTaskBase,
           URL = "task://localhost/dump_pread_ids" )
def dump_pread_ids(self):
    pread_db = fn( self.pread_db )
    pread_id_file = fn( self.pread_id_file )
    os.system("DBshow -n %s | tr -d '>' | awk '{print $1}' > %s" % (pread_db, pread_id_file) )

wf.addTask( dump_pread_ids )
wf.refreshTargets() # block

all_raw_las_files = {}
for las_fn in glob.glob( os.path.join( rawread_dir, "m*/raw_reads.*.las") ):
    idx = las_fn.split("/")[-1] # well, we will use regex someday to parse to get the number
    idx = int(idx.split(".")[1]) 
    las_file = makePypeLocalFile( las_fn )
    all_raw_las_files["r_las_%s" % idx] = las_file 

all_pread_las_files = {}
for las_fn in glob.glob( os.path.join( pread_dir, "m*/preads.*.las") ):
    idx = las_fn.split("/")[-1] # well, we will use regex someday to parse to get the number
    idx = int(idx.split(".")[1]) 
    las_file = makePypeLocalFile( las_fn )
    all_pread_las_files["p_las_%s" % idx] = las_file 



h_ctg_edges = makePypeLocalFile( os.path.join(asm_dir, "all_h_ctg_edges") )
p_ctg_edges = makePypeLocalFile( os.path.join(asm_dir, "all_p_ctg_edges") )
h_ctg_ids = makePypeLocalFile( os.path.join(asm_dir, "all_h_ctg_ids") )

inputs = { "rawread_id_file": rawread_id_file,
           "pread_id_file": pread_id_file,
           "h_ctg_edges": h_ctg_edges,
           "p_ctg_edges": p_ctg_edges,
           "h_ctg_ids": h_ctg_ids}

read_to_contig_map = makePypeLocalFile( os.path.join(read_map_dir, "read_to_contig_map") )

@PypeTask( inputs = inputs, 
           outputs = {"read_to_contig_map": read_to_contig_map}, 
           TaskType = PypeThreadTaskBase, 
           URL = "task://localhost/get_ctg_read_map" )

def generate_read_to_ctg_map(self):
    rawread_id_file = fn( self.rawread_id_file )
    pread_id_file = fn( self.pread_id_file )
    read_to_contig_map = fn( self.read_to_contig_map )
    
    pread_did_to_rid = open(pread_id_file).read().split("\n")
    rid_to_oid = open(rawread_id_file).read().split("\n")

    h_ctg_edges = fn( self.h_ctg_edges )
    p_ctg_edges = fn( self.p_ctg_edges )

    h_ctg_ids = set()
    with open(fn(self.h_ctg_ids)) as f:
        for row in f:
            row = row.strip()
            h_ctg_ids.add( row )

    pread_to_contigs = {}

    for fnanme in ( p_ctg_edges, h_ctg_edges):
        with open(fnanme) as f:
            for row in f:
                row = row.strip().split()
                ctg = row[0]
                if len(ctg.split("_")) > 1 and ctg not in h_ctg_ids:
                    continue
                n1 = row[1]
                n2 = row[2]
                pid1 = int(n1.split(":")[0])
                pid2 = int(n2.split(":")[0])
                rid1 = pread_did_to_rid[pid1].split("/")[1]
                rid2 = pread_did_to_rid[pid2].split("/")[1]
                rid1 = int(int(rid1)/10)
                rid2 = int(int(rid2)/10)
                oid1 = rid_to_oid[rid1]
                oid2 = rid_to_oid[rid2]
                k1 = (pid1, rid1, oid1)
                pread_to_contigs.setdefault( k1, set() )
                pread_to_contigs[ k1 ].add( ctg )
                k2 = (pid2, rid2, oid2)
                pread_to_contigs.setdefault( k2, set() )
                pread_to_contigs[ k2 ].add( ctg )

    with open(read_to_contig_map, "w") as f:
        for k in pread_to_contigs:
            pid, rid, oid = k
            for ctg in list(pread_to_contigs[ k ]):
                print >>f, "%09d %09d %s %s" % (pid, rid, oid, ctg)

wf.addTask( generate_read_to_ctg_map )

def dump_rawread_to_ctg(self):
    rawread_db = fn(self.rawread_db)
    rawread_id_file = fn(self.rawread_id_file)
    phased_read_file = fn(self.phased_reads)
    las_file = fn(self.las_file)
    rawread_to_contig_file = fn(self.rawread_to_contig_file)
    read_to_contig_map = fn(self.read_to_contig_map)
    rid_to_oid = open(rawread_id_file).read().split('\n')


    ovlp_data = []
    ovlp_count = 0
    longest_ovlp = 0
    a_id = None
    rid_to_contigs = {}
    
    with open(read_to_contig_map) as f:
        for row in f:
            row = row.strip().split()
            pid, rid, oid, ctg = row
            rid = int(rid)
            rid_to_contigs.setdefault( rid, (oid, set() ) )
            rid_to_contigs[ rid ][1].add( ctg )

    oid_to_phase = {}
    with open(phased_read_file) as f:
        for row in f:
            row = row.strip().split()
            ctg_id, block, phase = row[1:4]
            oid = row[6]
            block = int(block)
            phase = int(phase)
            oid_to_phase[ oid ] = (ctg_id, block, phase)


    with open(rawread_to_contig_file, "w") as f:
        ovlp_data = {}
        cur_read_id = None
        skip_rest = 0
        for row in sp.check_output(shlex.split("LA4Falcon -m %s %s " % (rawread_db, las_file)) ).splitlines():

            row = row.strip().split()
            t_id = int(row[1])
            q_id = int(row[0])
            if q_id != cur_read_id:
                if cur_read_id == None:
                    cur_read_id = q_id
                else:
                    if len(ovlp_data) == 0:
                        o_id = rid_to_oid[ cur_read_id ]
                        print >>f, "%09d %s %s %d %d %d %d" % (cur_read_id, o_id, "NA", 0, 0, 0, 0)
                    if len(ovlp_data) != 0:
                        ovlp_v = ovlp_data.values()
                        ovlp_v.sort()
                        rank = 0
                        for score, count, q_id_, o_id, ctg, in_ctg in ovlp_v:
                            print >> f, "%09d %s %s %d %d %d %d" % (q_id_, o_id, ctg, count, rank, score, in_ctg)
                            rank += 1
                    ovlp_data = {}
                    cur_read_id = q_id
                    skip_rest = 0

            if q_id in rid_to_contigs and len(ovlp_data) == 0: #if the query is already an edge of some contig....
                t_o_id, ctgs = rid_to_contigs[ q_id ]
                o_id = rid_to_oid[ q_id ]
                for ctg in list(ctgs):
                    ovlp_data.setdefault(ctg, [0, 0, q_id, o_id, ctg, 1])
                    ovlp_data[ctg][0] = -int(row[7]) 
                    ovlp_data[ctg][1] += 1
                    skip_rest = 1

            if skip_rest == 1:
                continue

            if t_id not in rid_to_contigs:
                continue

            q_phase = oid_to_phase.get( rid_to_oid[q_id], None )
            if q_phase != None:
                ctg_id, block, phase = q_phase
                if block != -1:
                    t_phase = oid_to_phase.get( rid_to_oid[t_id], None )
                    if t_phase != None:
                        if t_phase[0] == ctg_id and t_phase[1] == block and t_phase[2] != phase:
                            continue

            t_o_id, ctgs = rid_to_contigs[ t_id ]
            o_id = rid_to_oid[ q_id ]
            
            for ctg in list(ctgs):
                ovlp_data.setdefault(ctg, [0, 0, q_id, o_id, ctg, 0])
                ovlp_data[ctg][0] += int(row[2])
                ovlp_data[ctg][1] += 1

        if len(ovlp_data) != 0:
            ovlp_v = ovlp_data.values()
            ovlp_v.sort()
            rank = 0
            for score, count, q_id_, o_id, ctg, in_ctg in ovlp_v:
                print >> f, "%09d %s %s %d %d %d %d" % (q_id_, o_id, ctg, count, rank, score, in_ctg)
                rank += 1

def dump_pread_to_ctg(self):
    pread_db = fn( self.pread_db )
    rawread_id_file = fn( self.rawread_id_file )
    pread_id_file = fn( self.pread_id_file )
    phased_read_file = fn( self.phased_reads)
    read_to_contig_map = fn( self.read_to_contig_map )
    las_file = fn( self.las_file )
    pread_to_contig_file = fn( self.pread_to_contig_file )
    read_to_contig_map = fn( self.read_to_contig_map )
    
    pid_to_rid = open(pread_id_file).read().split("\n")
    rid_to_oid = open(rawread_id_file).read().split("\n")


    ovlp_data = []
    ovlp_count = 0
    longest_ovlp = 0
    a_id = None
    pid_to_contigs = {}
    
    with open(read_to_contig_map) as f:
        for row in f:
            row = row.strip().split()
            pid, rid, oid, ctg = row
            pid = int(pid)
            pid_to_contigs.setdefault( pid, (oid, set() ) )
            pid_to_contigs[ pid ][1].add( ctg )
            
    oid_to_phase = {}
    with open(phased_read_file) as f:
        for row in f:
            row = row.strip().split()
            ctg_id, block, phase = row[1:4]
            oid = row[6]
            block = int(block)
            phase = int(phase)
            oid_to_phase[ oid ] = (ctg_id, block, phase)

    with open(pread_to_contig_file, "w") as f:
        ovlp_data = {}
        cur_read_id = None
        skip_rest = 0
        for row in sp.check_output(shlex.split("LA4Falcon -mo %s %s " % (pread_db, las_file)) ).splitlines():

            row = row.strip().split()
            t_id = int(row[1])
            q_id = int(row[0])
            if q_id != cur_read_id:
                if cur_read_id == None:
                    cur_read_id = q_id
                else:
                    if len(ovlp_data) == 0:
                        rid = pid_to_rid[cur_read_id].split("/")[1]
                        rid = int(int(rid)/10)
                        o_id = rid_to_oid[ rid ]
                        print >>f, "%09d %s %s %d %d %d %d" % (cur_read_id, o_id, "NA", 0, 0, 0, 0)
                    else:
                        ovlp_v = ovlp_data.values()
                        ovlp_v.sort()
                        rank = 0
                        for score, count, q_id_, o_id, ctg, in_ctg in ovlp_v:
                            print >> f, "%09d %s %s %d %d %d %d" % (q_id_, o_id, ctg, count, rank, score, in_ctg)
                            rank += 1
                    ovlp_data = {}
                    cur_read_id = q_id
                    skip_rest = 0

            if q_id in pid_to_contigs and len(ovlp_data) == 0: #if the query is in some contig....
                t_o_id, ctgs = pid_to_contigs[ q_id ]
                rid = pid_to_rid[q_id].split("/")[1]
                rid = int(int(rid)/10)
                o_id = rid_to_oid[ rid ]
                for ctg in list(ctgs):
                    ovlp_data.setdefault(ctg, [0, 0, q_id, o_id, ctg, 1])
                    ovlp_data[ctg][0] = -int(row[7]) 
                    ovlp_data[ctg][1] += 1
                skip_rest = 1

            if skip_rest == 1:
                continue

            if t_id not in pid_to_contigs:
                continue
        
            q_rid = int( int(pid_to_rid[q_id].split("/")[1])/10 )
            q_phase = oid_to_phase.get( rid_to_oid[ q_rid ], None )
            
            if q_phase != None:
                ctg_id, block, phase = q_phase
                if block != -1:
                    t_rid = int( int(pid_to_rid[t_id].split("/")[1])/10 )
                    t_phase = oid_to_phase.get( rid_to_oid[ t_rid ], None )
                    if t_phase != None:
                        if t_phase[0] == ctg_id and t_phase[1] == block and t_phase[2] != phase:
                            continue

            t_o_id, ctgs = pid_to_contigs[ t_id ]
            rid = pid_to_rid[q_id].split("/")[1]
            rid = int(int(rid)/10)
            o_id = rid_to_oid[ rid ]
            
            for ctg in list(ctgs):
                ovlp_data.setdefault(ctg, [0, 0, q_id, o_id, ctg, 0])
                ovlp_data[ctg][0] += int(row[2])
                ovlp_data[ctg][1] += 1

        if len(ovlp_data) != 0:
            ovlp_v = ovlp_data.values()
            ovlp_v.sort()
            rank = 0
            for score, count, q_id_, o_id, ctg, in_ctg in ovlp_v:
                print >> f, "%09d %s %s %d %d %d %d" % (q_id_, o_id, ctg, count, rank, score, in_ctg)
                rank += 1


phased_reads =  makePypeLocalFile(os.path.join(asm_dir, "all_phased_reads"))


for las_key, las_file in all_raw_las_files.items():
    las_fn = fn(las_file)
    idx = las_fn.split("/")[-1] # well, we will use regex someday to parse to get the number
    idx = int(idx.split(".")[1]) 
    rawread_to_contig_file = makePypeLocalFile(os.path.join(read_map_dir, "rawread_to_contigs.%s" % idx))
    make_dump_rawread_to_ctg = PypeTask( inputs = { "las_file": las_file, 
                                                    "rawread_db": rawread_db, 
                                                    "read_to_contig_map": read_to_contig_map, 
                                                    "rawread_id_file": rawread_id_file,
                                                    "pread_id_file": pread_id_file,
                                                    "phased_reads" : phased_reads},
                                      outputs = { "rawread_to_contig_file": rawread_to_contig_file },
                                      TaskType = PypeThreadTaskBase,
                                      URL = "task://localhost/r_read_to_contigs.%s" % idx )
    dump_rawread_to_ctg_task = make_dump_rawread_to_ctg(dump_rawread_to_ctg)                           
    wf.addTask( dump_rawread_to_ctg_task )

for las_key, las_file in all_pread_las_files.items():
    las_fn = fn(las_file)
    idx = las_fn.split("/")[-1] # well, we will use regex someday to parse to get the number
    idx = int(idx.split(".")[1]) 
    pread_to_contig_file = makePypeLocalFile(os.path.join(read_map_dir, "pread_to_contigs.%s" % idx))
    make_dump_pread_to_ctg = PypeTask( inputs = { "las_file": las_file, 
                                                  "pread_db": pread_db, 
                                                  "read_to_contig_map": read_to_contig_map, 
                                                  "rawread_id_file": rawread_id_file,
                                                  "pread_id_file": pread_id_file,
                                                  "phased_reads": phased_reads},
                                      outputs = { "pread_to_contig_file": pread_to_contig_file },
                                      TaskType = PypeThreadTaskBase,
                                      URL = "task://localhost/pread_to_contigs.%s" % idx )
    dump_pread_to_ctg_task = make_dump_pread_to_ctg(dump_pread_to_ctg)                           
    wf.addTask( dump_pread_to_ctg_task )

wf.refreshTargets() # block
