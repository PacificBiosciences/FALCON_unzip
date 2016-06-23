from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader
import os
import networkx as nx
from multiprocessing import Pool
import argparse
import re
import sys

RCMAP = dict(zip("ACGTacgtNn-","TGCAtgcaNn-"))
## for shared memory usage
global p_asm_G
global h_asm_G
global all_rid_to_phase
global seqs

def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def reverse_end( node_id ):
    node_id, end = node_id.split(":")
    new_end = "B" if end == "E" else "E"
    return node_id + ":" + new_end


def load_sg_seq(all_read_ids, fasta_fn):

    seqs = {}
    # load all p-read name into memory
    f = FastaReader(fasta_fn)
    for r in f:
        if r.name not in all_read_ids:
            continue
        seqs[r.name] = r.sequence.upper()
    return seqs

def generate_haplotigs_for_ctg(input_):
   
    ctg_id, out_dir = input_
    global p_asm_G
    global h_asm_G
    global all_rid_to_phase
    global seqs
    arid_to_phase = all_rid_to_phase[ctg_id]

    mkdir( out_dir )

    ctg_G = p_asm_G.get_sg_for_ctg(ctg_id) 

    ctg_nodes = set(ctg_G.nodes())

    sg = nx.DiGraph()
    
    for v, w in ctg_G.edges():
        
        vrid = v[:9]
        wrid = w[:9]
            
        edge_data = p_asm_G.sg_edges[ (v, w) ]
        if edge_data[-1] != "G":
            continue

        vphase = arid_to_phase.get(vrid, (-1,0))
        wphase = arid_to_phase.get(wrid, (-1,0))
        if vphase[0] == wphase[0] and vphase[1] != wphase[1]:
            cross_phase = "Y"
        else:
            cross_phase = "N"

        sg.add_node( v, label= "%d_%d" % vphase, 
                        phase="%d_%d" % vphase, 
                        src="P" )

        sg.add_node( w, label= "%d_%d" % wphase, 
                        phase="%d_%d" % wphase, 
                        src="P" )

        sg.add_edge(v, w, src="OP", cross_phase = cross_phase)

        # we need to add the complimentary edges as the ctg_graph does not contain the dual edges
        rv = reverse_end(v)
        rw = reverse_end(w)
        sg.add_node( rv, label= "%d_%d" % vphase, 
                         phase="%d_%d" % vphase, 
                         src="P" )

        sg.add_node( rw, label= "%d_%d" % wphase, 
                         phase="%d_%d" % wphase, 
                         src="P" )

        sg.add_edge(rw, rv, src="OP", cross_phase = cross_phase)

    PG_nodes = set(sg.nodes())
    PG_edges = set(sg.edges())

    for v, w in h_asm_G.sg_edges:
        
        vrid = v[:9]
        wrid = w[:9]

        if vrid not in arid_to_phase:
            continue
        if wrid not in arid_to_phase:
            continue
        
        if (v, w) in PG_edges:
            if p_asm_G.sg_edges[(v,w)][-1] == "G":
                continue

        edge_data = h_asm_G.sg_edges[ (v, w) ]

        if edge_data[-1] != "G":
            continue

        cross_phase = "N"
        if v not in PG_nodes:
            sg.add_node( v, label= "%d_%d" % arid_to_phase[vrid], 
                            phase="%d_%d" % arid_to_phase[vrid], 
                            src="H" )

        if w not in PG_nodes:
            sg.add_node( w, label= "%d_%d" % arid_to_phase[wrid], 
                            phase="%d_%d" % arid_to_phase[wrid], 
                            src="H" )

        sg.add_edge(v, w, src="H", cross_phase = cross_phase)

        rv = reverse_end(v)
        rw = reverse_end(w)
        if rv not in PG_nodes:
            sg.add_node( rv, label= "%d_%d" % arid_to_phase[vrid], 
                             phase="%d_%d" % arid_to_phase[vrid], 
                             src="H" )

        if rw not in PG_nodes:
            sg.add_node( rw, label= "%d_%d" % arid_to_phase[wrid], 
                             phase="%d_%d" % arid_to_phase[wrid], 
                             src="H" )

        sg.add_edge(rw, rv, src="H", cross_phase = cross_phase)

    sg0 = sg.copy()
    for v, w in h_asm_G.sg_edges:
        vrid = v[:9]
        wrid = w[:9]
        if vrid not in arid_to_phase:
            continue
        if wrid not in arid_to_phase:
            continue
        
        if (v, w) in PG_edges:
            if p_asm_G.sg_edges[(v,w)][-1] == "G":
                continue

        edge_data = h_asm_G.sg_edges[ (v, w) ]

        if sg0.in_degree(w) == 0:
            cross_phase = "Y"
            if v not in PG_nodes:
                sg.add_node( v, label= "%d_%d" % arid_to_phase[vrid], 
                                phase="%d_%d" % arid_to_phase[vrid], 
                                src="H" )

            if w not in PG_nodes:
                sg.add_node( w, label= "%d_%d" % arid_to_phase[wrid], 
                                phase="%d_%d" % arid_to_phase[wrid], 
                                src="H" )

            sg.add_edge(v, w, src="ext", cross_phase = cross_phase)

            rv = reverse_end(v)
            rw = reverse_end(w)
            if rv not in PG_nodes:
                sg.add_node( rv, label= "%d_%d" % arid_to_phase[vrid], 
                                 phase="%d_%d" % arid_to_phase[vrid], 
                                 src="H" )

            if rw not in PG_nodes:
                sg.add_node( rw, label= "%d_%d" % arid_to_phase[wrid], 
                                 phase="%d_%d" % arid_to_phase[wrid], 
                                 src="H" )

            sg.add_edge(rw, rv, src="ext", cross_phase = cross_phase)

        if sg0.out_degree(v) == 0:
            cross_phase = "Y"
            if v not in PG_nodes:
                sg.add_node( v, label= "%d_%d" % arid_to_phase[vrid], 
                                phase="%d_%d" % arid_to_phase[vrid], 
                                src="H" )

            if w not in PG_nodes:
                sg.add_node( w, label= "%d_%d" % arid_to_phase[wrid], 
                                phase="%d_%d" % arid_to_phase[wrid], 
                                src="H" )

            sg.add_edge(v, w, src="ext", cross_phase = cross_phase)

            rv = reverse_end(v)
            rw = reverse_end(w)
            if rv not in PG_nodes:
                sg.add_node( rv, label= "%d_%d" % arid_to_phase[vrid], 
                                 phase="%d_%d" % arid_to_phase[vrid], 
                                 src="H" )

            if rw not in PG_nodes:
                sg.add_node( rw, label= "%d_%d" % arid_to_phase[wrid], 
                                 phase="%d_%d" % arid_to_phase[wrid], 
                                 src="H" )

            sg.add_edge(rw, rv, src="ext", cross_phase = cross_phase)

    sg2 = sg.copy()
    ctg_nodes_r = set([ reverse_end(v) for v in list(ctg_nodes) ])
    for v, w in ctg_G.edges():
        sg2.remove_edge(v, w)
        rv, rw = reverse_end(v), reverse_end(w)
        sg2.remove_edge(rw, rv)
    for v in sg2.nodes():
        if sg2.out_degree(v) == 0 and sg2.in_degree(v) == 0:
            sg2.remove_node(v)

    nodes_to_remove = set()
    edges_to_remove = set()
    for sub_g in nx.weakly_connected_component_subgraphs(sg2):
        sub_g_nodes = set(sub_g.nodes())
        if len(sub_g_nodes & ctg_nodes_r) > 0 and len(sub_g_nodes & ctg_nodes) > 0:
            # remove cross edge
            sources = [n for n in sub_g.nodes() if sub_g.in_degree(n) == 0 or n in ctg_nodes or n in ctg_nodes_r ]
            sinks = [n for n in sub_g.nodes() if sub_g.out_degree(n) == 0 or n in ctg_nodes or n in ctg_nodes_r ]
            edges_to_keep = set()
            for v in sources:
                for w in sinks:
                    path = []
                    if v in ctg_nodes and w not in ctg_nodes_r:
                        try:
                            path = nx.shortest_path( sub_g, v, w ) 
                        except nx.exception.NetworkXNoPath:
                            path = []
                    elif v not in ctg_nodes and w in ctg_nodes_r:
                        try:
                            path = nx.shortest_path( sub_g, v, w )
                        except nx.exception.NetworkXNoPath:
                            path = []

                    if len(path) >= 2:
                        v1 = path[0]
                        for w1 in path[1:]:
                            edges_to_keep.add( (v1, w1) )
                            rv1, rw1 = reverse_end(v1), reverse_end(w1)
                            edges_to_keep.add( (rw1, rv1) )
                            v1 = w1
            for v, w in sub_g.edges():
                if (v, w) not in edges_to_keep:
                    edges_to_remove.add( (v, w) )
                    rv, rw = reverse_end(v), reverse_end(w)
                    edges_to_remove.add( (rw, rv) )


        if len(sub_g_nodes & ctg_nodes_r) == 0 and len(sub_g_nodes & ctg_nodes) == 0:
            nodes_to_remove.update( sub_g_nodes )
            nodes_to_remove.update( set( [reverse_end(v) for v in list(sub_g_nodes)] ) )

    for v, w in list(edges_to_remove):
        sg.remove_edge(v, w)

    for v in nodes_to_remove:
        sg.remove_node(v)

    for v in sg.nodes():
        if sg.out_degree(v) == 0 and sg.in_degree(v) == 0:
            sg.remove_node(v)

    #nx.write_gexf(sg, "full_g.gexf")
    
    s_node = p_asm_G.ctg_data[ctg_id][5][0][0]
    t_node = p_asm_G.ctg_data[ctg_id][5][-1][-1]

    for v, w in sg.edges():
        phase0 = sg.node[v]["phase"].split("_")
        phase1 = sg.node[w]["phase"].split("_")
        if phase0 == phase1:
            sg[v][w]["weight"] = 10
            sg[v][w]["score"] = 1
            sg[v][w]["label"] = "type0" 
        else:
            if phase0[0] == phase1[0]:
                sg[v][w]["weight"] = 1
                sg[v][w]["score"] = 100000
                sg[v][w]["label"] = "type1"
            else:
                sg[v][w]["weight"] = 5
                sg[v][w]["score"] = 50
                sg[v][w]["label"] = "type2"


    sg2 = sg.copy()
    edge_to_remove = set()
    for v, w in sg2.edges():
        if sg2[v][w]["src"] == "ext":
            edge_to_remove.add( (v, w) )
            rv, rw = reverse_end(v), reverse_end(w)
            edge_to_remove.add( (rw, rv) )

        if sg2.node[v]["phase"] ==  sg2.node[w]["phase"]:
            continue
        flag1 = 0
        flag2 = 0
        for e in sg2.out_edges(v):
            if sg2.node[e[0]]["phase"] ==  sg2.node[e[1]]["phase"]:
                flag1 = 1
                break
        if flag1 == 1:
            for e in sg2.in_edges(w):
                if sg2.node[e[0]]["phase"] ==  sg2.node[e[1]]["phase"]:
                    flag2 = 1
                    break
        if flag2 == 1:
            edge_to_remove.add( (v, w) )
            rv, rw = reverse_end(v), reverse_end(w)
            edge_to_remove.add( (rw, rv) )


    for v, w in list(edge_to_remove):
        sg2.remove_edge(v, w)
    try: 
        s_path = nx.shortest_path(sg2, source=s_node, target=t_node, weight="score")
    except nx.exception.NetworkXNoPath:
        s_path = nx.shortest_path(sg, source=s_node, target=t_node, weight="score")

    s_path_edges = [] 
    for i in xrange(len(s_path)-1):
        v = s_path[i]
        w = s_path[i+1]
        sg[v][w]["weight"] = 15
        s_path_edges.append( (v,w) )

    s_path_edge_set = set(s_path_edges)


    
    #output the updated primary contig
    p_tig_path = open(os.path.join(out_dir, "p_ctg_path.%s" % ctg_id),"w")
    p_tig_fa = open(os.path.join(out_dir, "p_ctg.%s.fa" % ctg_id),"w")
    edges_to_remove1 = set()
    edges_to_remove2 = set()
    with open(os.path.join(out_dir, "p_ctg_edges.%s" % ctg_id), "w") as f:
        seq = []
        for v, w in s_path_edges:
            sg[v][w]["h_edge"] = 1
            vrid = v.split(":")[0]
            wrid = w.split(":")[0]
            vphase = arid_to_phase.get(vrid, (-1,0))
            wphase = arid_to_phase.get(wrid, (-1,0))
            print >>f, "%s" % ctg_id, v, w, sg[v][w]["cross_phase"], sg[v][w]["src"], vphase[0], vphase[1], wphase[0], wphase[1]

            if sg.edge[v][w]["src"] == "OP":
                edge_data = p_asm_G.sg_edges[ (v,w) ]
            else:
                edge_data = h_asm_G.sg_edges[ (v,w) ]

            seq_id, s, t = edge_data[0]
            if s < t:
                seq.append(seqs[ seq_id ][ s:t ])
            else:
                seq.append("".join([ RCMAP[c] for c in seqs[ seq_id ][ s:t:-1 ] ]))
            print >>p_tig_path, "%s" % ctg_id, v, w, seq_id, s, t, edge_data[1], edge_data[2], "%d %d" % arid_to_phase.get(seq_id, (-1,0))
            sg[v][w]["tig_id"] = "%s" % ctg_id

            rv, rw = reverse_end(v), reverse_end(w)
            edges_to_remove1.add( (v, w) )
            edges_to_remove2.add( (rw, rv) )

        print >> p_tig_fa, ">%s" % ctg_id
        print >> p_tig_fa, "".join(seq)

    p_tig_fa.close()
    p_tig_path.close()



    sg2 = sg.copy()
    reachable1 = nx.descendants(sg2, s_node)
    sg2_r = sg2.reverse()
    reachable2 = nx.descendants(sg2_r, t_node)

    reachable_all = reachable1 | reachable2
    reachable_both = reachable1 & reachable2


    for v, w in list(edges_to_remove2 | edges_to_remove1):
        sg2.remove_edge( v, w )

    for v, w in sg2.edges():
        if sg2[v][w]["cross_phase"] == "Y":
            sg2.remove_edge( v, w )

    for v in sg2.nodes():
        if v not in reachable_all:
            sg2.remove_node(v)

    for v in sg2.nodes():
        if sg2.out_degree(v) == 0 and sg2.in_degree(v) == 0:
            sg2.remove_node(v)
            continue
        if v in reachable_both:
            sg2.node[v]["reachable"] = 1
        else:
            sg2.node[v]["reachable"] = 0
        
    dump_graph = False # the code segement below is useful for showing the graph
    if dump_graph == True:
        nx.write_gexf(sg2, "%s_1.gexf" % ctg_id)
    
    p_path_nodes = set(s_path)
    p_path_rc_nodes = set( [reverse_end(v) for v in s_path] )

    sg2_nodes = set(sg2.nodes())
    for v in p_asm_G.get_sg_for_ctg(ctg_id).nodes():
        rv = reverse_end(v)
        p_path_rc_nodes.add( rv )
        if rv in sg2_nodes:
            sg2.remove_node(rv)

    
    h_tig_path = open(os.path.join(out_dir, "h_ctg_path.%s" % ctg_id),"w")
    h_tig_fa = open(os.path.join(out_dir, "h_ctg_all.%s.fa" % ctg_id),"w")
    edges_to_remove = set()

    labelled_node = set()
    with open(os.path.join(out_dir, "h_ctg_edges.%s" % ctg_id),"w") as f:
        h_tig_id = 1
        h_paths = {}
        #print "number of components:", len([tmp for tmp in nx.weakly_connected_component_subgraphs(sg2)])
        for sub_hg_0 in nx.weakly_connected_component_subgraphs(sg2):
            sub_hg = sub_hg_0.copy()
            while sub_hg.size() > 5:
                #print "sub_hg size:", len(sub_hg.nodes())
                sources = [n for n in sub_hg.nodes() if sub_hg.in_degree(n) != 1 ]
                sinks = [n for n in sub_hg.nodes() if sub_hg.out_degree(n) != 1 ]
                

                #print "number of sources", len(sources),  sources
                #print "number of sinks", len(sinks), sinks
                if len(sources) == 0 and len(sinks) == 0: #TODO, the rest of the sub-graph are circles, we need to break and print warnning message
                    break

                longest = [] 

                eliminated_sinks = set()
                s_longest = {}
                for s in sources:
                    #print "test source",s, len(eliminated_sinks)
                    if s in labelled_node:
                        continue
                    s_path = []
                    for t in sinks:
                        if t in eliminated_sinks:
                            continue
                        try:
                            path = nx.shortest_path(sub_hg, s, t, weight="score")
                            #print "test path len:", len(path), s, t
                        except nx.exception.NetworkXNoPath:
                            path = []
                            continue
                        s_path.append( [ path, t ] )
                    s_path.sort(key = lambda x: -len(x[0]))
                    if len(s_path) == 0:
                        continue
                    s_longest[s] = s_path[0][0]
                    if len(s_longest[s]) > len(longest):
                        longest = s_longest[s]
                        #print "s longest", longest[0], longest[-1], len(longest)
                    for path, t in s_path[1:]:
                        eliminated_sinks.add(t)
                        #print "elimated t", t
                            

                if len(longest) == 0:
                    break

                s = longest[0]
                t = longest[-1]
                h_paths[ ( s, t ) ] = longest
                
                labelled_node.add(s)
                rs = reverse_end(s)
                labelled_node.add(rs)

                for v in longest:
                    sub_hg.remove_node(v)

        for s, t in h_paths:
            longest = h_paths[ (s, t) ]
            #print "number of node in path", s,t,len(longest) 
            seq = []
            for v, w in zip(longest[:-1], longest[1:]):
                sg[v][w]["h_edge"] = 1
                if sg.edge[v][w]["src"] == "OP":
                    edge_data = p_asm_G.sg_edges[ (v,w) ]
                else:
                    edge_data = h_asm_G.sg_edges[ (v,w) ]
                vrid = v.split(":")[0]
                wrid = w.split(":")[0]
                vphase = arid_to_phase.get(vrid, (-1,0))
                wphase = arid_to_phase.get(wrid, (-1,0))
                print >>f, "%s_%03d" % (ctg_id, h_tig_id), v, w, sg[v][w]["cross_phase"], sg[v][w]["src"], vphase[0], vphase[1], wphase[0], wphase[1]

                if sg.edge[v][w]["src"] == "OP":
                    edge_data = p_asm_G.sg_edges[ (v,w) ]
                else:
                    edge_data = h_asm_G.sg_edges[ (v,w) ]

                seq_id, sp, tp = edge_data[0]
                if sp < tp:
                    seq.append(seqs[ seq_id ][ sp:tp ])
                else:
                    seq.append("".join([ RCMAP[c] for c in seqs[ seq_id ][ sp:tp:-1 ] ]))
                print >> h_tig_path, "%s_%03d" % (ctg_id, h_tig_id), v, w, seq_id, sp, tp, edge_data[1], edge_data[2], "%d %d" % arid_to_phase.get(seq_id, (-1,0))
                sg[v][w]["tig_id"] = "%s_%03d" % (ctg_id, h_tig_id)

                rv, rw = reverse_end(v), reverse_end(w)
                edges_to_remove.add( (v, w) )
                edges_to_remove.add( (rw, rv) )

            print >> h_tig_fa, ">%s_%03d" % (ctg_id, h_tig_id)
            print >> h_tig_fa, "".join(seq)
            h_tig_id += 1


    h_tig_fa.close()
    h_tig_path.close()

    dump_graph = False  # the code segement below is useful for showing the graph
    if dump_graph == True:
        for v, w in sg.edges():
            if "h_edge" not in sg[v][w]:
                sg[v][w]["h_edge"] = 0
            if v in reachable_all:
                sg.node[v]["reachable"] = 1
            else:
                sg.node[v]["reachable"] = 0
            if w in reachable_all:
                sg.node[w]["reachable"] = 1
            else:
                sg.node[w]["reachable"] = 0

        nx.write_gexf(sg, "%s_0.gexf" % ctg_id)



def parse_args(argv):

    parser = argparse.ArgumentParser(description='layout haplotigs from primary assembly graph and phased aseembly graph')

    parser.add_argument('--fc_asm_path', type=str, help='path to the primary Falcon assembly output directory', required=True)
    parser.add_argument('--fc_hasm_path', type=str, help='path to the phased Falcon assembly output directory', required=True)
    parser.add_argument('--ctg_id', type=str, help='contig identifier in the bam file', default = "all", required=True)
    parser.add_argument('--base_dir', type=str, default="./", help='the output base_dir, default to current working directory')
    parser.add_argument('--rid_phase_map', type=str, help="path to the file that encode the relationship of the read id to phase blocks", required=True)
    parser.add_argument('--fasta', type=str, help="sequence file of the p-reads", required=True)
    args = parser.parse_args(argv[1:])

    return args

def main(argv=sys.argv):

    # make life easier for now. will refactor it out if possible
    global all_rid_to_phase
    global p_asm_G
    global h_asm_G
    global all_rid_to_phase
    global seqs

    args = parse_args(argv)
    fc_asm_path = args.fc_asm_path
    fc_hasm_path = args.fc_hasm_path
    ctg_id = args.ctg_id
    base_dir = args.base_dir
    fasta_fn = args.fasta
    
    p_asm_G = AsmGraph(os.path.join(fc_asm_path, "sg_edges_list"), 
                       os.path.join(fc_asm_path, "utg_data"),
                       os.path.join(fc_asm_path, "ctg_paths") )

    h_asm_G = AsmGraph( os.path.join(fc_hasm_path, "sg_edges_list"), 
                        os.path.join(fc_hasm_path, "utg_data"), 
                        os.path.join(fc_hasm_path, "ctg_paths") )


    all_rid_to_phase = {}

    all_read_ids = set()
    with open(args.rid_phase_map) as f:
        for row in f:
            row = row.strip().split()
            all_rid_to_phase.setdefault( row[1], {} )
            all_rid_to_phase[row[1]][row[0]] = (int(row[2]), int(row[3]))
            all_read_ids.add(row[0])

    for v, w in p_asm_G.sg_edges:
        if p_asm_G.sg_edges[ (v, w) ][-1] != "G":
            continue
        v = v.split(":")[0]
        w = w.split(":")[0]
        all_read_ids.add(v)
        all_read_ids.add(w)

    for v, w in h_asm_G.sg_edges:
        if h_asm_G.sg_edges[ (v, w) ][-1] != "G":
            continue
        v = v.split(":")[0]
        w = w.split(":")[0]
        all_read_ids.add(v)
        all_read_ids.add(w)

    seqs = load_sg_seq(all_read_ids, fasta_fn)

    if ctg_id == "all":
        ctg_id_list = p_asm_G.ctg_data.keys()
    else:
        ctg_id_list = [ctg_id]

    exe_list = []
    for ctg_id in ctg_id_list:
        if ctg_id[-1] != "F":
            continue
        if ctg_id not in all_rid_to_phase:
            continue
        exe_list.append( (ctg_id, os.path.join(".", ctg_id)) )

    exec_pool = Pool(4)  #TODO, make this configurable
    exec_pool.map( generate_haplotigs_for_ctg, exe_list)
    #map( generate_haplotigs_for_ctg, exe_list)
