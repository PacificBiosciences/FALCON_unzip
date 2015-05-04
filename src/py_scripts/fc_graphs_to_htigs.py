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
    phase_connection = {}
    for v, w in ctg_G.edges():
        
        vrid = v[:9]
        wrid = w[:9]
            
        edge_data = p_asm_G.sg_edges[ (v, w) ]
        if edge_data[-1] != "G":
            continue
        if vrid not in arid_to_phase:
            continue
        if wrid not in arid_to_phase:
            continue
        sg.add_node( v, label= "%d_%d" % arid_to_phase[vrid], phase="%d_%d" % arid_to_phase[vrid], src="P" )
        sg.add_node( w, label= "%d_%d" % arid_to_phase[wrid], phase="%d_%d" % arid_to_phase[wrid], src="P" )
        if arid_to_phase[vrid][0] == arid_to_phase[wrid][0] and arid_to_phase[vrid][1] != arid_to_phase[wrid][1]:
            cross_phase = "Y"
        else:
            cross_phase = "N"
        if arid_to_phase[vrid][0] != arid_to_phase[wrid][0]:
            phase_connection.setdefault( arid_to_phase[vrid], set() )
            if arid_to_phase[vrid][0] != -1 and arid_to_phase[wrid][0] != -1:
                phase_connection[ arid_to_phase[vrid] ].add( arid_to_phase[wrid] )

        sg.add_edge(v, w, src="P", cross_phase = cross_phase)

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

        sg.add_edge(v, w, src="H", cross_phase = "N")
        

    s_node = p_asm_G.ctg_data[ctg_id][5][0][0]
    h_ctg_G = sg.subgraph(nx.descendants(sg, s_node))

    h_ctg_G_0 = h_ctg_G.copy()

    h_ctg_G_orig = h_ctg_G.copy()
    for e in h_ctg_G_orig.edges():
        if h_ctg_G_orig.edge[e[0]][e[1]]['cross_phase'] == "Y":
                h_ctg_G.remove_edge(e[0], e[1])


    h_ctg_G_orig = h_ctg_G.copy()
    for n in h_ctg_G_orig.nodes():
        for v, w in h_ctg_G_orig.out_edges(n):
            vrid = v.split(":")[0]
            wrid = w.split(":")[0]
            if arid_to_phase[vrid] == arid_to_phase[wrid]: #if the edge connects the same phase, keep the edge
                continue

            phase = arid_to_phase[vrid]
            anti_phase = list(phase)
            anti_phase[1] = 1 - anti_phase[1]
            anti_phase = tuple(anti_phase)

            #keep the edge if the connection to the next phase is consistent
            if len(phase_connection.get(phase,set()) ) == 1 and len(phase_connection.get(anti_phase,set())) == 1:
                p1 = list(phase_connection[ phase  ])[0]
                p2 = list(phase_connection[ anti_phase ])[0]
                if p1[0] == p2[0] and p1[1] != p2[1]:
                    continue
            h_ctg_G.remove_edge(v, w) #remove all other edges with different phase 
            
    for n in h_ctg_G.nodes():
        if h_ctg_G.in_degree(n) == 0 and h_ctg_G.out_degree(n) == 0:
            h_ctg_G.remove_node(n)

    nx.write_gexf(h_ctg_G, "%s_1.gexf" % ctg_id)

    h_edges = h_ctg_G.edges()
    for e in h_ctg_G_0.edges():
        if e in h_edges:    
            h_ctg_G_0.edge[e[0]][e[1]]["h_edge"] = "Y"
        else:
            h_ctg_G_0.edge[e[0]][e[1]]["h_edge"] = "N"

    nx.write_gexf(h_ctg_G_0, "%s_0.gexf" % ctg_id)


    all_read_ids = set()
    for v, w in h_ctg_G.edges():
        v = v.split(":")[0]
        w = w.split(":")[0]
        all_read_ids.add(v)
        all_read_ids.add(w)

    seqs = load_sg_seq(all_read_ids, fasta_fn)


    h_tig_fa = open( os.path.join(base_dir, "h_ctg.fa"), "w")
    h_tig_path = open( os.path.join(base_dir, "h_ctg_path"), "w")
    h_tig_edge = open( os.path.join(base_dir, "h_ctg_edge"), "w")

    h_tig_id = 0
    for sub_hg in nx.weakly_connected_component_subgraphs(h_ctg_G):
        sources = [n for n in sub_hg.nodes() if sub_hg.in_degree(n) == 0]
        sinks = [n for n in sub_hg.nodes() if sub_hg.out_degree(n) == 0]
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
        if len(longest) < 5:
            continue
        seq = []
        for v, w in zip(longest[:-1], longest[1:]):
            if h_ctg_G.edge[v][w]["src"] == "P":
                edge_data = p_asm_G.sg_edges[ (v,w) ]
            else:
                edge_data = h_asm_G.sg_edges[ (v,w) ]

            seq_id, s, t = edge_data[0]
            if s < t:
                seq.append(seqs[ seq_id ][ s:t ])
            else:
                seq.append("".join([ RCMAP[c] for c in seqs[ seq_id ][ s:t:-1 ] ]))
            print >>h_tig_path, "%s_H%03d" % (ctg_id, h_tig_id), v, w, seq_id, s, t, edge_data[1], edge_data[2], "%d %d" % arid_to_phase[seq_id] 
        print >> h_tig_fa, ">%s_H%03d" % (ctg_id, h_tig_id)
        print >> h_tig_fa, "".join(seq)
        h_tig_id += 1
        
        for v, w in sub_hg.edges():
            if h_ctg_G.edge[v][w]["src"] == "P":
                edge_data = p_asm_G.sg_edges[ (v,w) ]
            else:
                edge_data = h_asm_G.sg_edges[ (v,w) ]

            seq_id, s, t = edge_data[0]
            print >>h_tig_edge, "%s_H%03d" % (ctg_id, h_tig_id), v, w, seq_id, s, t, edge_data[1], edge_data[2], "%d %d" % arid_to_phase[seq_id]


    h_tig_fa.close()
    h_tig_edge.close()
    h_tig_path.close()

