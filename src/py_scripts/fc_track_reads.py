from pypeflow.common import * 
from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
from pypeflow.task import PypeTask, PypeThreadTaskBase, PypeTaskBase
from pypeflow.controller import PypeWorkflow, PypeMPWorkflow, PypeThreadWorkflow
from falcon_kit.FastaReader import FastaReader
from falcon_kit.fc_asm_graph import AsmGraph
import glob
import sys
import subprocess as sp
import shlex
import os

rawread_dir = os.path.abspath( "./0-rawreads" )
pread_dir = os.path.abspath( "./1-preads_ovl" )
asm_dir = os.path.abspath( "./2-asm-falcon" )

PypeMPWorkflow.setNumThreadAllowed(12, 12)
wf = PypeMPWorkflow()

def dump_rawread_map(self):
    rawread_db = fn( self.rawread_db )
    las_file = fn( self.las_file )
    ovlp_map_id_file = fn( self.ovlp_map_id_file )

    ovlp_data = []
    ovlp_count = 0
    longest_ovlp = 0
    a_id = None
    with open(ovlp_map_id_file, "w") as f:
        for row in sp.check_output(shlex.split("LA4Falcon -mo %s %s " % (rawread_db, las_file)) ).splitlines():
            row = row.strip().split()
            if row[-1] == "contained":
                continue

            if row[0] != a_id:
                if a_id != None and len(ovlp_data) > 0:
                    ovlp_data.sort()
                    out_ids = []
                    for d in ovlp_data[:100]:
                        out_ids.append( d[1][1] ) #query id
                    print >>f, a_id, ovlp_count, " ".join(out_ids)
                ovlp_count = 0
                longest_ovlp = 0
                ovlp_data = []
                a_id = row[0]

            if int(row[9]) > 500 or int(row[11]) - int(row[10]) > 500: #only output contained relationship
                continue
            ovlp_len  = int(row[6]) - int(row[5])
            ovlp_data.append( (-ovlp_len, row) )
            ovlp_count += 1

        # last element    
        if a_id != None and len(ovlp_data) > 0:
            ovlp_data.sort()
            out_ids = []
            for d in ovlp_data[:100]:
                out_ids.append( d[1][1] ) #query id
            print >>f, row[0], ovlp_count, " ".join(out_ids)

def dump_pread_map(self):
    pread_db = fn( self.pread_db )
    las_file = fn( self.las_file )
    ovlp_map_id_file = fn( self.ovlp_map_id_file )
    a_id = None
    with open(ovlp_map_id_file, "w") as f:
        for row in sp.check_output(shlex.split("LA4Falcon -mo %s %s " % (pread_db, las_file)) ).splitlines():
            row = row.strip().split()

            if row[0] != a_id:
                if a_id != None and len(ovlp_data) > 0:
                    #ovlp_data.sort()
                    out_ids = []
                    for d in ovlp_data:
                        out_ids.append( d[1][1] ) #query id
                    print >>f, a_id, ovlp_count, " ".join(out_ids)
                ovlp_count = 0
                longest_ovlp = 0
                ovlp_data = []
                a_id = row[0]

            if int(row[9]) > 5 or int(row[11]) - int(row[10]) > 5: #only output contained relationship
                continue
            ovlp_len  = int(row[6]) - int(row[5])
            ovlp_data.append( (-ovlp_len, row) )
            ovlp_count += 1

        # last element    
        if a_id != None and len(ovlp_data) > 0:
            ovlp_data.sort()
            out_ids = []
            for d in ovlp_data:
                out_ids.append( d[1][1] ) #query id
            print >>f, row[0], ovlp_count, " ".join(out_ids)


rawread_db = makePypeLocalFile( os.path.join( rawread_dir, "raw_reads.db" ) )
rawread_id_file = makePypeLocalFile( os.path.join( rawread_dir, "raw_reads_ids" ) )

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
pread_id_file = makePypeLocalFile( os.path.join( pread_dir, "preads_ids" ) )

@PypeTask( inputs = {"pread_db": pread_db}, 
           outputs =  {"pread_id_file": pread_id_file},
           TaskType = PypeThreadTaskBase,
           URL = "task://localhost/dump_pread_ids" )
def dump_pread_ids(self):
    pread_db = fn( self.pread_db )
    pread_id_file = fn( self.pread_id_file )
    os.system("DBshow -n %s | tr -d '>' | awk '{print $1}' > %s" % (pread_db, pread_id_file) )

wf.addTask( dump_pread_ids )

all_raw_las_files = {}
all_r_ovlp_map_files = {}
for las_fn in glob.glob( os.path.join( rawread_dir, "raw_reads.*.las") ):
    idx = las_fn.split("/")[-1] # well, we will use regex someday to parse to get the number
    idx = int(idx.split(".")[1]) 

    las_file = makePypeLocalFile( las_fn )
    all_raw_las_files["r_las_%s" % idx] = las_file 
    ovlp_map_id_file  = makePypeLocalFile( os.path.join( rawread_dir, "ovlp_map.%s" % idx ) ) 
    all_r_ovlp_map_files["r_idmap_%s" % idx ] = ovlp_map_id_file
    make_rawread_map_task = PypeTask( inputs = { "las_file": las_file, "rawread_db": rawread_db },
                                      outputs = { "ovlp_map_id_file": ovlp_map_id_file },
                                      TaskType = PypeThreadTaskBase,
                                      URL = "task://localhost/r_ovlp_map.%s" % idx )
    rawread_map_task = make_rawread_map_task(dump_rawread_map)                            
    wf.addTask( rawread_map_task )

all_p_las_files = {}
all_p_ovlp_map_files = {}
for las_fn in glob.glob( os.path.join( pread_dir, "preads.*.las") ):
    idx = las_fn.split("/")[-1] # well, we will use regex someday to parse to get the number
    idx = int(idx.split(".")[1]) 

    las_file = makePypeLocalFile( las_fn )
    all_p_las_files["p_las_%s" % idx] = las_file 
    ovlp_map_id_file  = makePypeLocalFile( os.path.join( pread_dir, "ovlp_map.%s" % idx ) ) 
    all_p_ovlp_map_files["p_idmap_%s" % idx ] = ovlp_map_id_file
    make_pread_map_task = PypeTask( inputs = { "las_file": las_file, "pread_db": pread_db },
                                    outputs = { "ovlp_map_id_file": ovlp_map_id_file },
                                    TaskType = PypeThreadTaskBase,
                                    URL = "task://localhost/p_ovlp_map.%s" % idx )
    pread_map_task = make_pread_map_task(dump_pread_map)                            
    wf.addTask( pread_map_task )

wf.refreshTargets() # block

# need new workflow
PypeMPWorkflow.setNumThreadAllowed(1, 1)
wf = PypeMPWorkflow()

p_ctg_tiling_path = makePypeLocalFile( os.path.join(asm_dir, "p_ctg_tiling_path") )
a_ctg_tiling_path = makePypeLocalFile( os.path.join(asm_dir, "a_ctg_tiling_path") )
sg_edges_list = makePypeLocalFile( os.path.join(asm_dir, "sg_edges_list") )
utg_data = makePypeLocalFile( os.path.join(asm_dir, "utg_data") )
ctg_paths = makePypeLocalFile( os.path.join(asm_dir, "ctg_paths") )

inputs = { "rawread_id_file": rawread_id_file,
           "pread_id_file": pread_id_file,
           "sg_edges_list": sg_edges_list,
           "utg_data": utg_data,
           "ctg_paths": ctg_paths }

parameters = { "p_ovlp_files": [ fn(fobj) for fobj in all_p_ovlp_map_files.values() ],
               "r_ovlp_files": [ fn(fobj) for fobj in all_r_ovlp_map_files.values() ] }  # put these files as parameters to avoid tracking them

contig_to_read_map = makePypeLocalFile( os.path.join(asm_dir, "contig_to_read_map") )

@PypeTask( inputs = inputs, 
           outputs = {"contig_to_read_map": contig_to_read_map}, 
           parameters = parameters,  
           TaskType = PypeThreadTaskBase, 
           URL = "task://localhost/get_ctg_read_map" )
def gen_ctg_to_read_map(self):
    rawread_id_file = fn( self.rawread_id_file )
    pread_id_file = fn( self.pread_id_file )
    contig_to_read_map = fn( self.contig_to_read_map )
    p_ovlp_map_files = self.parameters["p_ovlp_files"]
    r_ovlp_map_files = self.parameters["r_ovlp_files"]
    p_ovlp_map_files.sort()
    r_ovlp_map_files.sort()
    tiling_files = [ fn(self.p_ctg_tiling_path), fn(self.a_ctg_tiling_path) ]
    
    pread_did_to_rid = open(pread_id_file).read().split("\n")
    rid_to_oid = open(rawread_id_file).read().split("\n")
    
    
    overlap_pread = {}
    for ovlp_fn in p_ovlp_map_files:
        with open(ovlp_fn) as f:
            for row in f:
                row = row.strip().split()
                overlap_pread[int(row[0])] = set([ int(c) for c in row[2:] ])

    overlap_rawread = {}
    for ovlp_fn in r_ovlp_map_files:
        with open(ovlp_fn) as f:
            for row in f:
                row = row.strip().split()
                overlap_rawread[int(row[0])] = set([ int(c) for c in row[2:] ])

    ctg_to_preads = {}
    rid_set = set()
    asm_G = AsmGraph(fn(self.sg_edges_list), 
                       fn(self.utg_data),
                       fn(self.ctg_paths) )


    with open(contig_to_read_map, "w") as f:
        for ctg in asm_G.ctg_data:
            rid_set = set()
            ctg_to_preads = {}
            if ctg[-1] == "R":
                continue
            ctg_g = asm_G.get_sg_for_ctg(ctg)
            for n in ctg_g.nodes():
                frg0 = int(n.split(":")[0])

                ctg_to_preads.setdefault( ctg, set() )

                rid = pread_did_to_rid[frg0].split("/")[1]
                rid = int(rid[:-1])
                oid = rid_to_oid[rid]
                ctg_to_preads[ctg].add((0, frg0, rid, oid))
                rid_set.add(frg0)

                for frg in list(overlap_pread.get(frg0, set())):
                    rid = pread_did_to_rid[frg].split("/")[1]
                    rid = int(int(rid) / 10)
                    if rid in rid_set:
                        continue
                    oid = rid_to_oid[rid]
                    ctg_to_preads[ctg].add((1, frg, rid, oid))
                    rid_set.add(rid)


            for i in range(2,4):
                for r in list(ctg_to_preads[ctg]):
                    class_, pid, rid0, oid = r

                    for rid in list(overlap_rawread.get(rid0,set())):
                        if rid in rid_set:
                            continue
                        oid = rid_to_oid[rid]
                        ctg_to_preads[ctg].add((i ,rid0, rid, oid))
                        rid_set.add(rid)


            for ctg in ctg_to_preads:
                for r in list(ctg_to_preads[ctg]):
                    print >>f, ctg, r[0], "%09d" % r[1], "%09d" % r[2], r[3]

wf.addTask( gen_ctg_to_read_map )
wf.refreshTargets()
