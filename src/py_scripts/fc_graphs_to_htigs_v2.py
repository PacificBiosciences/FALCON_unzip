from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader
import os
import networkx as nx

RCMAP = dict(zip("ACGTacgtNn-","TGCAtgcaNn-"))


def load_sg_seq(all_read_ids, fasta_fn):

    seqs = {}
    # load all p-read name into memory
    f = FastaReader(fasta_fn)
    for r in f:
        if r.name not in all_read_ids:
            continue
        seqs[r.name] = r.sequence.upper()
    return seqs

if __name__ == "__main__":
    import argparse
    import re
    parser = argparse.ArgumentParser(description='layout haplotigs from primary assembly graph and phased aseembly graph')
    parser.add_argument('--fc_asm_path', type=str, help='path to the primary Falcon assembly output directory', required=True)
    parser.add_argument('--fc_hasm_path', type=str, help='path to the phased Falcon assembly output directory', required=True)
    parser.add_argument('--ctg_id', type=str, help='contig identifier in the bam file', required=True)
    parser.add_argument('--base_dir', type=str, default="./", help='the output base_dir, default to current working directory')
    parser.add_argument('--rid_phase_map', type=str, help="path to the file that encode the relationship of the read id to phase blocks", required=True)
    parser.add_argument('--fasta', type=str, help="sequence file of the p-reads", required=True)

    args = parser.parse_args()
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

    ctg_G = p_asm_G.get_sg_for_ctg(ctg_id) 

    ctg_nodes = set(ctg_G.nodes())

    arid_to_phase = {}

    with open(args.rid_phase_map) as f:
        for row in f:
            row = row.strip().split()
            arid_to_phase[row[0]] = (int(row[1]), int(row[2]))

    sg = nx.DiGraph()
    
    for v, w in ctg_G.edges():
        
        vrid = v[:9]
        wrid = w[:9]
            
        edge_data = p_asm_G.sg_edges[ (v, w) ]
        if edge_data[-1] != "G":
            continue
        """
        if vrid not in arid_to_phase:
            continue
        if wrid not in arid_to_phase:
            continue
        """
        vphase = arid_to_phase.get(vrid, (-1,0))
        wphase = arid_to_phase.get(wrid, (-1,0))
        sg.add_node( v, label= "%d_%d" % vphase, phase="%d_%d" % vphase, src="P" )
        sg.add_node( w, label= "%d_%d" % wphase, phase="%d_%d" % wphase, src="P" )
        if vphase[0] == wphase[0] and vphase[1] != wphase[1]:
            cross_phase = "Y"
        else:
            cross_phase = "N"

        sg.add_edge(v, w, src="OP", cross_phase = cross_phase)

    PG_nodes = set(sg.nodes())
    PG_edges = set(sg.edges())
    

    for v, w in h_asm_G.sg_edges:
        
        vrid = v[:9]
        wrid = w[:9]

        if (v, w) in PG_edges:
            if p_asm_G.sg_edges[(v,w)][-1] == "G":
                continue

        edge_data = h_asm_G.sg_edges[ (v, w) ]

        if edge_data[-1] != "G":
            continue

        if v not in PG_nodes:
            sg.add_node( v, label= "%d_%d" % arid_to_phase[vrid], phase="%d_%d" % arid_to_phase[vrid], src="H" )
        if w not in PG_nodes:
            sg.add_node( w, label= "%d_%d" % arid_to_phase[wrid], phase="%d_%d" % arid_to_phase[wrid], src="H" )
        cross_phase = "N"

        sg.add_edge(v, w, src="H", cross_phase = cross_phase)
        

    s_node = p_asm_G.ctg_data[ctg_id][5][0][0]
    t_node = p_asm_G.ctg_data[ctg_id][5][-1][-1]
    for g in nx.weakly_connected_component_subgraphs(sg):
        if s_node in g.nodes():
            h_ctg_G = g
            break

    sg = h_ctg_G

    for v, w in sg.edges():
        phase0 = sg.node[v]["phase"].split("_")
        phase1 = sg.node[w]["phase"].split("_")
        if phase0 == phase1:
            sg[v][w]["weight"] = 10
            sg[v][w]["score"] = 1
            sg[v][w]["label"] = "t0"
        else:
            if phase0[0] == phase1[0]:
                sg[v][w]["weight"] = 1
                sg[v][w]["score"] = 100000
                sg[v][w]["label"] = "t1"
            else:
                sg[v][w]["weight"] = 5
                sg[v][w]["score"] = 50
                sg[v][w]["label"] = "t2"

    node_to_remove = set()
    for v in sg.nodes():
        try:
            nx.shortest_path(sg, source=s_node, target=v)
        except nx.NetworkXNoPath:
            node_to_remove.add(v)
        try:
            nx.shortest_path(sg, source=v, target=t_node)
        except nx.NetworkXNoPath:
            node_to_remove.add(v)
    for v in list(node_to_remove):
        sg.remove_node(v)



    sg2 = sg.copy()
    edge_to_remove = set()
    for v, w in sg2.edges():
        if sg2.node[v]["phase"] ==  sg2.node[w]["phase"]:
            continue
        flag1 = 0
        flag2 = 0
        for e in sg2.out_edges(v):
            if sg2.node[e[0]]["phase"] ==  sg2.node[e[1]]["phase"]:
                flag1 = 1
        if flag1 == 1:
            for e in sg2.in_edges(w):
                if sg2.node[e[0]]["phase"] ==  sg2.node[e[1]]["phase"]:
                    flag2 = 1
        if flag2 == 1:
            edge_to_remove.add( (v, w) )


    for v, w in list(edge_to_remove):
        sg2.remove_edge(v,w)
        sg[v][w]["remove"] = 1
    
    s_path = nx.shortest_path(sg2, source=s_node, target=t_node, weight="score")

    s_path_edges = [] 
    for i in xrange(len(s_path)-1):
        v = s_path[i]
        w = s_path[i+1]
        sg[v][w]["weight"] = 15
        print v,w,sg[v][w]["score"]
        s_path_edges.append( (v,w) )

    s_path_edge_set = set(s_path_edges)
    for v, w in s_path_edges:
        sg2.remove_edge(v,w)
        
    for v in sg2.nodes():
        if sg2.out_degree(v) == 0 and sg2.in_degree(v) == 0:
            sg2.remove_node(v)

    for v, w in sg.edges():
        if "remove" not in sg[v][w]:
            sg[v][w]["remove"] = 0



    all_read_ids = set()
    for v, w in sg.edges():
        v = v.split(":")[0]
        w = w.split(":")[0]
        all_read_ids.add(v)
        all_read_ids.add(w)

    seqs = load_sg_seq(all_read_ids, fasta_fn)
    #output the updated primary contig

    h_tig_path = open("h_ctg_path","w")
    h_tig_fa = open("h_ctg.fa","w")
    p_path_nodes = set(s_path)
    with open("h_ctg_edges","w") as f:
        h_tig_id = 1
        for sub_hg in nx.weakly_connected_component_subgraphs(sg2):
            sources = [n for n in sub_hg.nodes() if sub_hg.in_degree(n) == 0 and n in p_path_nodes]
            sinks = [n for n in sub_hg.nodes() if sub_hg.out_degree(n) == 0 and n in p_path_nodes]
            longest = [] 
            for s in sources:
                for t in sinks:
                    try:
                        path = nx.shortest_path(sub_hg, s, t)
                    except nx.exception.NetworkXNoPath:
                        path = []
                        pass
                    if len(path) > len(longest):
                        longest = path

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

                seq_id, s, t = edge_data[0]
                if s < t:
                    seq.append(seqs[ seq_id ][ s:t ])
                else:
                    seq.append("".join([ RCMAP[c] for c in seqs[ seq_id ][ s:t:-1 ] ]))
                print >> h_tig_path, "%s_%03d" % (ctg_id, h_tig_id), v, w, seq_id, s, t, edge_data[1], edge_data[2], "%d %d" % arid_to_phase.get(seq_id, (-1,0))
            print >> h_tig_fa, ">%s_%04d" % (ctg_id, h_tig_id)
            print >> h_tig_fa, "".join(seq)
            h_tig_id += 1

    p_tig_path = open("p_ctg_path","w")
    p_tig_fa = open("p_ctg.fa","w")
    with open("p_ctg_edges","w") as f:
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
        print >> p_tig_fa, ">%s" % ctg_id
        print >> p_tig_fa, "".join(seq)



    for v, w in sg.edges():
        if "h_edge" not in sg[v][w]:
            sg[v][w]["h_edge"] = 0

    nx.write_gexf(sg, "%s_0.gexf" % ctg_id)
