from pypeflow.common import * 
from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
from pypeflow.task import PypeTask, PypeThreadTaskBase, PypeTaskBase
from pypeflow.controller import PypeWorkflow, PypeThreadWorkflow
from falcon_kit.FastaReader import FastaReader
import subprocess, shlex
import os
import re
import argparse
import sys

cigar_re = r"(\d+)([MIDNSHP=X])"

def make_het_call(self):
    
    bam_fn = fn(self.bam_file)
    ctg_id = self.parameters["ctg_id"]
    ref_seq = self.parameters["ref_seq"]
    base_dir = self.parameters["base_dir"]
    vmap_fn = fn(self.vmap_file)
    vpos_fn = fn(self.vpos_file)
    q_id_map_fn = fn(self.q_id_map_file)

    p = subprocess.Popen(shlex.split("samtools view %s %s" % (bam_fn, ctg_id) ), stdout=subprocess.PIPE)
    pileup = {}
    q_id_map = {}
    q_max_id = 0
    q_id = 0
    q_name_to_id = {}

    try:
        os.makedirs("%s/%s" % (base_dir, ctg_id))
    except OSError:
        pass

    vmap = open(vmap_fn, "w")
    vpos = open(vpos_fn, "w")

    for l in p.stdout:
        l = l.strip().split()
        if l[0][0] == "@":
            continue

        QNAME = l[0]
        if QNAME not in q_name_to_id:
            q_id = q_max_id
            q_name_to_id[QNAME] = q_id
            q_max_id += 1

        q_id = q_name_to_id[QNAME]
        q_id_map[q_id] = QNAME
        FLAG = int(l[1])
        RNAME = l[2]
        POS = int(l[3]) - 1 # convert to zero base
        CIGAR = l[5]
        SEQ = l[9]
        rp = POS
        qp = 0

        skip_base = 0
        total_aln_pos = 0
        for m in re.finditer(cigar_re, CIGAR):
            adv = int(m.group(1))
            total_aln_pos += adv

            if m.group(2)  == "S":
                skip_base += adv

        if 1.0 - 1.0 * skip_base / total_aln_pos < 0.1:
            continue
        if total_aln_pos < 2000:
            continue

        for m in re.finditer(cigar_re, CIGAR):
            adv = int(m.group(1))
            if m.group(2) == "S":
                qp += adv
            if m.group(2) == "M":
                matches = []
                for i in range(adv):
                    matches.append( (rp, SEQ[qp]) )
                    rp += 1
                    qp += 1
                matches = matches[1:-1]
                for pos, b in  matches:
                    pileup.setdefault(pos, {})
                    pileup[pos].setdefault(b, [])
                    pileup[pos][b].append(q_id)
            elif m.group(2) == "I":
                for i in range(adv):
                    qp += 1
            elif m.group(2) == "D":
                for i in range(adv):
                    rp += 1

        pos_k = pileup.keys()
        pos_k.sort()
        th = 0.25
        for pos in pos_k:
            if pos < POS:
                if len(pileup[pos]) < 2:
                    del pileup[pos]
                    continue
                base_count = [] 
                total_count = 0
                for b in ["A", "C", "G", "T"]:
                    count = len(pileup[pos].get(b,[]))
                    base_count.append( (count, b) )
                    total_count += count
                if total_count < 10:
                    del pileup[pos]
                    continue

                base_count.sort()
                base_count.reverse()
                p0 = 1.0 *  base_count[0][0] / total_count
                p1 = 1.0 *  base_count[1][0] / total_count
                if p0 < 1.0 - th and p1 > th:  
                    b0 = base_count[0][1]
                    b1 = base_count[1][1]
                    ref_base = ref_seq[pos]
                    print >> vpos, pos+1, ref_base, total_count, " ".join(["%s %d" % (x[1], x[0]) for x in base_count])
                    for q_id_ in pileup[pos][b0]:
                        print >> vmap, pos+1, ref_base, b0, q_id_
                    for q_id_ in pileup[pos][b1]:
                        print >> vmap, pos+1, ref_base, b1, q_id_ 
                del pileup[pos]


    q_id_map_f = open(q_id_map_fn, "w")
    for q_id, q_name in q_id_map.items():
        print >> q_id_map_f, q_id, q_name


def generate_association_table(self):

    vmap_fn = fn(self.vmap_file)
    atable_fn = fn(self.atable_file)
    ctg_id = self.parameters["ctg_id"]
    base_dir = self.parameters["base_dir"]

    vmap = {}
    v_positions = []
    
    with open(vmap_fn) as f:
        for l in f:
            l = l.strip().split()
            pos = int(l[0])
            ref_b = l[1]
            v_b = l[2]
            q_id = int(l[3])
            if (pos, ref_b) not in vmap:
                v_positions.append( (pos, ref_b) )
            vmap.setdefault( (pos, ref_b), {} )
            vmap[ (pos, ref_b) ].setdefault(v_b, [])
            vmap[ (pos, ref_b) ][v_b].append( q_id )


    #xary = []
    #yary = []
    with open(atable_fn, "w") as out_f:
        for i1 in xrange(len(v_positions)):
            link_count = 0
            for i2 in xrange(i1+1, len(v_positions)):
                pos1, rb1 = v_positions[i1]
                pos2, rb2 = v_positions[i2]
                if pos2 - pos1 > (1 << 16):
                    continue
                ct = {}
                p1table = []
                p2table = []
                s1 = 0
                list1 = vmap[ (pos1, rb1) ].items()
                for b1, qids1 in list1:
                    p1table.append( (b1, len(qids1) ) )
                    s1 += len(qids1)
                 
                s2 = 0
                list2 = vmap[ (pos2, rb2) ].items()
                for b2, qids2 in list2:
                    p2table.append( (b2, len(qids2) ) )
                    s2 += len(qids2)
                
                total_s = 0
                for b1, qids1 in list1:
                    for b2, qids2 in list2:
                        s = len(set(qids1) & set(qids2))
                        ct[(b1,b2)] = s
                        total_s += s
                if total_s < 6:
                    continue
                
                b11 = p1table[0][0]
                b12 = p1table[1][0]
                b21 = p2table[0][0]
                b22 = p2table[1][0]
                print >> out_f, pos1, b11, b12, pos2, b21, b22, ct[(b11,b21)], ct[(b11,b22)], ct[(b12,b21)], ct[(b12,b22)] 
            

                #xary.append(pos1)
                #yary.append(pos2)
                link_count += 1
                if link_count > 500:
                    break

def get_score( c_score, pos1, pos2, s1, s2 ):
    if pos1 > pos2:
        pos1, pos2 = pos2, pos1
        s1, s2 = s2, s1
    b11, b12 = s1
    b21, b22 = s2
    return c_score[ (pos1, pos2) ][ (b11+b21, b12+b22) ]

def get_phased_blocks(self):
    vmap_fn = fn(self.vmap_file)
    atable_fn = fn(self.atable_file)
    p_variant_fn = fn(self.phased_variant_file)

    left_connect = {}
    right_connect = {}

    c_score = {}
    states = {}
    positions = set() 



    ref_base = {}
    with open(vmap_fn) as f:
        for l in f:
            l = l.strip().split()
            pos = int(l[0])
            ref_b = l[1]
            v_b = l[2]
            q_id = int(l[3])
            ref_base[pos] = ref_b

    with open(atable_fn) as f:
        for l in f:
            l = l.strip().split()
            pos1, b11, b12, pos2, b21, b22, s11, s12, s21, s22 = l
            s11, s12, s21, s22 = int(s11), int(s12), int(s21), int(s22)
            if abs(s11+s22-s12-s21) < 6:
                continue
            pos1 = int(pos1)
            pos2 = int(pos2)
            positions.add(pos1)
            positions.add(pos2)
            right_connect.setdefault(pos1, [])
            right_connect[pos1].append(pos2)
            left_connect.setdefault(pos2, [])
            left_connect[pos2].append(pos1)
            c_score[ (pos1, pos2) ] = { (b11+b21, b12+b22): s11 + s22, (b12+b22, b11+b21): s11 + s22,
                                        (b12+b21, b11+b22): s12 + s21, (b11+b22, b12+b21): s12 + s21 }


            if pos1 not in states:
                st1 = (b11, b12)
                st2 = (b12, b11)
                score1 = 0
                score2 = 0
                for pp in left_connect.get(pos1,[]):
                    if pp in states:
                        st0 = states[pp]
                    else:
                        continue
                    score1 += get_score( c_score, pp, pos1, st0, st1 )
                    score2 += get_score( c_score, pp, pos1, st0, st2 )

                for pp in right_connect.get(pos1,[]):
                    if pp in states:
                        st0 = states[pp]
                    else:
                        continue
                    score1 += get_score( c_score, pos1, pp, st1, st0 )
                    score2 += get_score( c_score, pos1, pp, st2, st0 )

                if score1 >= score2:
                    states[pos1] = st1
                else:
                    states[pos1] = st2

            if pos2 not in states:
                st1 = (b21, b22)
                st2 = (b22, b21)
                score1 = 0
                score2 = 0
                for pp in left_connect.get(pos2,[]):
                    if pp in states:
                        st0 = states[pp]
                    else:
                        continue
                    score1 += get_score( c_score, pp, pos2, st0, st1 )
                    score2 += get_score( c_score, pp, pos2, st0, st2 )

                for pp in right_connect.get(pos2,[]):
                    if pp in states:
                        st0 = states[pp]
                    else:
                        continue
                    score1 += get_score( c_score, pos2, pp, st1, st0 )
                    score2 += get_score( c_score, pos2, pp, st2, st0 )

                if score1 >= score2:
                    states[pos2] = st1
                else:
                    states[pos2] = st2

    positions = list(positions)
    positions.sort()


    iter_count = 0
    while 1:
        iter_count += 1
        if iter_count > 10:
            break
        update_count = 0
        for p in positions:
            b1, b2 = states[p]
            st1 = (b1, b2)
            st2 = (b2, b1)

            score1 = 0
            score2 = 0
            for pp in left_connect.get(p,[]):
                st0 = states[pp]
                score1 += get_score( c_score, pp, p, st0 ,st1)
                score2 += get_score( c_score, pp, p, st0, st2)

            #for pp in right_connect.get(p,[]):
            #    st0 = states[pp]
            #    score1 += get_score( c_score, p, pp, st1 ,st0)
            #    score2 += get_score( c_score, p, pp, st2, st0)

            if score1 >= score2:
                states[p] = st1
            else:
                states[p] = st2
                update_count += 1
        if update_count == 0:
            break


    right_extent = {}
    right_score = {}
    left_extent = {}
    left_score = {}

    
    for p in positions:

        left_extent[p] = p
        left_score[p] = 0
        if p in left_connect:
            left = p
            st0 = states[p]
            st0_ = st0[1], st0[0]
            for pp in left_connect[p]:
                st1 = states[pp]
                s = get_score( c_score, pp, p, st1, st0)
                s_ = get_score( c_score, pp, p, st1, st0_) 
                left_score[p] += s - s_
                if s - s_ > 0 and pp < left:
                    left = pp
            left_extent[p] = left
                
        right_extent[p] = p
        right_score[p] = 0
        if p in right_connect:
            right = p
            st0 = states[p]
            st0_ = st0[1], st0[0]
            for pp in right_connect[p]:
                st1 = states[pp]
                s = get_score( c_score, p, pp, st0, st1)
                s_ = get_score( c_score, p, pp, st0_, st1) 
                right_score[p] += s - s_
                if s - s_ > 0 and pp > right:
                    right = pp
            right_extent[p] = right
        
        


    phase_block_id = 1
    phase_blocks = {}
    pb = []

    max_right_ext = 0
    for p in positions:
        if right_score[p] < 10 or left_score[p] < 10:
            continue
        b1, b2 = states[p]
        if max_right_ext < left_extent[p]:
            if len(pb) > 3:
                phase_blocks[phase_block_id] = pb
                phase_block_id += 1
            pb = []
        pb.append( (p, b1, b2) )
        if right_extent[p] > max_right_ext:
            max_right_ext =  right_extent[p]
    if len(pb) > 3:
        phase_blocks[phase_block_id] = pb
    else:
        phase_block_id -= 1

    
    with open(p_variant_fn, "w") as out_f:
        for pid in xrange(1, phase_block_id+1):
            if len(phase_blocks[pid]) == 0:
                continue
            min_ = min( [x[0] for x in phase_blocks[pid]] )
            max_ = max( [x[0] for x in phase_blocks[pid]] )
            
            print >>out_f, "P", pid, min_, max_, max_ - min_, len(phase_blocks[pid]), 1.0 * (max_-min_)/len(phase_blocks[pid]) 
            for p, b1, b2 in phase_blocks[pid]:
                rb = ref_base[p]
                print >>out_f, "V", pid, p, "%d_%s_%s" % (p,rb,b1), "%d_%s_%s" % (p,rb,b2), left_extent[p], right_extent[p], left_score[p], right_score[p] 

def get_phased_reads(self):

    q_id_map_fn = fn(self.q_id_map_file)
    vmap_fn = fn(self.vmap_file)
    p_variant_fn = fn(self.phased_variant_file)
    parameters = self.parameters

    ctg_id = parameters["ctg_id"]

    phased_read_fn = fn(self.phased_read_file) 

    rid_map = {}
    with open(q_id_map_fn) as f:
        for l in f:
            l = l.strip().split()
            rid_map[int(l[0])] = l[1]


    read_to_variants = {}
    variant_to_reads = {}
    with open(vmap_fn) as f: 
        for l in f:
            l = l.strip().split()
            variant = "_".join(l[:3])
            read_id = int(l[3])
            read_to_variants.setdefault(read_id, set())
            read_to_variants[read_id].add(variant)
            variant_to_reads.setdefault(variant, set())
            variant_to_reads[variant].add(read_id)


    variant_to_phase = {}
    with open(p_variant_fn) as f:
        for l in f:
            """line format example: V 1 6854 6854_A_A 6854_A_G 6854 22781"""
            l = l.strip().split()
            if l[0] != "V":
                continue
            pb_id = int(l[1])
            variant_to_phase[ l[3] ] = (pb_id, 0)
            variant_to_phase[ l[4] ] = (pb_id, 1)
    
    with open(phased_read_fn, "w") as out_f:
        for r in read_to_variants:
            vl = {}
            pl = set()
            for v in list( read_to_variants[r] ):
                if v in variant_to_phase:
                    p = variant_to_phase[v]
                    vl[ p ] = vl.get(p, 0) + 1
                    pl.add(p[0])
            pl = list(pl)
            pl.sort()
            for p in pl:
                if vl.get( (p,0), 0) - vl.get( (p,1), 0) > 1:
                    print >> out_f, r, ctg_id, p, 0, vl.get( (p,0), 0), vl.get( (p,1), 0), rid_map[r]
                elif vl.get( (p,1), 0) - vl.get( (p,0), 0) > 1:
                    print >> out_f, r, ctg_id, p, 1, vl.get( (p,0), 0), vl.get( (p,1), 0), rid_map[r]

def phasing(args):
    bam_fn = args.bam
    fasta_fn = args.fasta
    ctg_id = args.ctg_id
    base_dir = args.base_dir
    
    ref_seq = "" 
    for r in FastaReader(fasta_fn):
        rid = r.name.split()[0]
        if rid != ctg_id:
            continue
        ref_seq = r.sequence.upper()

    PypeThreadWorkflow.setNumThreadAllowed(1, 1)
    wf = PypeThreadWorkflow()



    bam_file = makePypeLocalFile(bam_fn)
    vmap_file = makePypeLocalFile( os.path.join(base_dir, ctg_id, "variant_map") )
    vpos_file = makePypeLocalFile( os.path.join(base_dir, ctg_id, "variant_pos") )
    q_id_map_file = makePypeLocalFile( os.path.join(base_dir, ctg_id, "q_id_map") )
    parameters = {}
    parameters["ctg_id"] = ctg_id
    parameters["ref_seq"] = ref_seq
    parameters["base_dir"] = base_dir
    
    make_het_call_task = PypeTask( inputs = { "bam_file": bam_file },
                         outputs = { "vmap_file": vmap_file, "vpos_file": vpos_file, "q_id_map_file": q_id_map_file },
                         parameters = parameters,
                         TaskType = PypeThreadTaskBase,
                         URL = "task://localhost/het_call") (make_het_call)

    wf.addTasks([make_het_call_task])




    atable_file = makePypeLocalFile( os.path.join(base_dir, ctg_id, "atable") )
    parameters = {}
    parameters["ctg_id"] = ctg_id
    parameters["base_dir"] = base_dir
    generate_association_table_task = PypeTask( inputs = { "vmap_file": vmap_file },
                                      outputs = { "atable_file": atable_file },
                                      parameters = parameters,
                                      TaskType = PypeThreadTaskBase,
                                      URL = "task://localhost/g_atable") (generate_association_table)

    wf.addTasks([generate_association_table_task])




    phased_variant_file = makePypeLocalFile( os.path.join(base_dir, ctg_id, "phased_variants") )
    get_phased_blocks_task = PypeTask( inputs = { "vmap_file": vmap_file, "atable_file": atable_file },
                                      outputs = { "phased_variant_file": phased_variant_file },
                                      TaskType = PypeThreadTaskBase,
                                      URL = "task://localhost/get_phased_blocks") (get_phased_blocks)
    wf.addTasks([get_phased_blocks_task])




    phased_read_file = makePypeLocalFile( os.path.join(base_dir, ctg_id, "phased_reads") )
    get_phased_reads_task = PypeTask( inputs = { "vmap_file": vmap_file, 
                                                 "q_id_map_file": q_id_map_file, 
                                                 "phased_variant_file": phased_variant_file },
                                      outputs = { "phased_read_file": phased_read_file },
                                      parameters = {"ctg_id": ctg_id},
                                      TaskType = PypeThreadTaskBase,
                                      URL = "task://localhost/get_phased_reads") (get_phased_reads)
    wf.addTasks([get_phased_reads_task])
    

    wf.refreshTargets() 
    #with open("fc_phasing_wf.dot", "w") as f:
    #    print >>f, wf.graphvizDot 

def parse_args(argv):
    parser = argparse.ArgumentParser(description='phasing variants and reads from a bam file')
    # we can run this in parallel mode in the furture
    #parser.add_argument('--n_core', type=int, default=4,
    #                    help='number of processes used for generating consensus')
    parser.add_argument('--bam', type=str, help='path to sorted bam file', required=True)
    parser.add_argument('--fasta', type=str, help='path to the fasta file of contain the contig', required=True)
    parser.add_argument('--ctg_id', type=str, help='contig identifier in the bam file', required=True)
    parser.add_argument('--base_dir', type=str, default="./", help='the output base_dir, default to current working directory')


    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    phasing(args)

