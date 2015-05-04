from falcon_kit.fc_asm_graph import AsmGraph
from falcon_kit.FastaReader import FastaReader
import networkx as nx

RCMAP = dict(zip("ACGTacgtNn-","TGCAtgcaNn-"))

CVI_set0 = set()
with open("../data/CVI_bax_input.fofn") as f:
    for row in f:
        CVI_set0.add(row.strip().split("/")[-1].split(".")[0])

COL_set0 = set()
with open("../data/COL_bax_input.fofn") as f:
    for row in f:
        COL_set0.add(row.strip().split("/")[-1].split(".")[0])

p_asm_G = AsmGraph("../2-asm-falcon/sg_edges_list", "../2-asm-falcon/utg_data", "../2-asm-falcon/ctg_paths")
h_asm_G = AsmGraph("../5-diploid-graph/sg_edges_list", "../5-diploid-graph/utg_data", "../5-diploid-graph/ctg_paths")

ctg = "000000F"
ctg_G = p_asm_G.get_sg_for_ctg(ctg) 

ctg_nodes = set(ctg_G.nodes())
rid_to_phase = {}
with open("../4-contig-phasing/phased_reads") as f:
    for row in f:
        row = row.strip().split()
        rid_to_phase[row[5]] = (int(row[1]), int(row[2]))

arid_to_phase = {}
CVI_set = set()
COL_set = set()
with open("../2-asm-falcon/contig_to_read_map") as f:
    for row in f:
        row = row.strip().split()
        ctg_id = row[0]
        if not ctg_id.startswith("000000F"):
            continue
        if int(row[1]) > 1: #not preads
            continue
        phase = rid_to_phase.get( row[4], (-1, 0) )
        arid_to_phase["%09d" % int(row[2])] = phase

        if row[4].split("/")[0] in CVI_set0:
            CVI_set.add( "%09d" % int(row[2]) )
        if row[4].split("/")[0] in COL_set0:
            COL_set.add( "%09d" % int(row[2]) )

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
    if vrid in CVI_set:
        v_cvi = "CVI"
    else:
        assert vrid in COL_set
        v_cvi = "COL"
    if wrid in CVI_set:
        w_cvi = "CVI"
    else:
        assert wrid in COL_set
        w_cvi = "COL"
    sg.add_node( v, label= "%d_%d" % arid_to_phase[vrid], phase="%d_%d" % arid_to_phase[vrid], src="P", CVI = v_cvi )
    sg.add_node( w, label= "%d_%d" % arid_to_phase[wrid], phase="%d_%d" % arid_to_phase[wrid], src="P", CVI = w_cvi )
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

    if vrid in CVI_set:
        v_cvi = "CVI"
    else:
        assert vrid in COL_set
        v_cvi = "COL"
    if wrid in CVI_set:
        w_cvi = "CVI"
    else:
        assert wrid in COL_set
        w_cvi = "COL"

    if v not in PG_nodes:
        sg.add_node( v, label= "%d_%d" % arid_to_phase[vrid], phase="%d_%d" % arid_to_phase[vrid], src="H", CVI = v_cvi )
    if w not in PG_nodes:
        sg.add_node( w, label= "%d_%d" % arid_to_phase[wrid], phase="%d_%d" % arid_to_phase[wrid], src="H", CVI = w_cvi )

    sg.add_edge(v, w, src="H", cross_phase = "N")
    

s_node = p_asm_G.ctg_data["000000F"][5][0][0]
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

nx.write_gexf(h_ctg_G,"000000F_1.gexf")

h_edges = h_ctg_G.edges()
for e in h_ctg_G_0.edges():
    if e in h_edges:    
        h_ctg_G_0.edge[e[0]][e[1]]["hedge"] = "Y"
    else:
        h_ctg_G_0.edge[e[0]][e[1]]["hedge"] = "N"


nx.write_gexf(h_ctg_G_0,"000000F.gexf")

def load_sg_seq(all_read_ids):

    fasta_fn = "./preads4falcon.fasta"


    seqs = {}
    # load all p-read name into memory
    f = FastaReader(fasta_fn)
    for r in f:
        if r.name not in all_read_ids:
            continue
        seqs[r.name] = r.sequence.upper()
    return seqs

all_read_ids = set()
for v, w in h_ctg_G.edges():
    v = v.split(":")[0]
    w = w.split(":")[0]
    all_read_ids.add(v)
    all_read_ids.add(w)

seqs = load_sg_seq(all_read_ids)


h_tig_fa = open("h_ctg.fa", "w")
h_tig_path = open("h_ctg_path","w")
h_tig_edge = open("h_ctg_edge","w")

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
        print >>h_tig_path, "000000F_H%03d" % h_tig_id, v, w, seq_id, s, t, edge_data[1], edge_data[2], "%d %d" % arid_to_phase[seq_id] 
    print >> h_tig_fa, ">000000F_H%03d" % h_tig_id
    print >> h_tig_fa, "".join(seq)
    h_tig_id += 1
    
    for v, w in sub_hg.edges():
        if h_ctg_G.edge[v][w]["src"] == "P":
            edge_data = p_asm_G.sg_edges[ (v,w) ]
        else:
            edge_data = h_asm_G.sg_edges[ (v,w) ]

        seq_id, s, t = edge_data[0]
        print >>h_tig_edge, "000000F_H%03d" % h_tig_id, v, w, seq_id, s, t, edge_data[1], edge_data[2], "%d %d" % arid_to_phase[seq_id]


h_tig_fa.close()
h_tig_edge.close()
h_tig_path.close()

