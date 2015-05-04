import subprocess, shlex
import re
from falcon_kit.FastaReader import FastaReader


cigar_re = r"(\d+)([MIDNSHP=X])"

ref_seq = {}
for r in FastaReader("./000000F_ref.fa"):
    ref_seq[r.name.split()[0]] = r.sequence.upper()


p = subprocess.Popen(shlex.split("samtools_ view ../3-contig_mapping/aln_sorted.bam 000000F" ), stdout=subprocess.PIPE)
pileup = {}
q_id_map = {}
q_id = 0

vmap2 = open("vmap2","w")
vpos2 = open("vpos2","w")

for l in p.stdout:
    l = l.strip().split()
    if l[0][0] == "@":
        continue

    QNAME = l[0]
    q_id_map[q_id] = QNAME
    FLAG = int(l[1])
    RNAME = l[2]
    POS = int(l[3]) - 1 # convert to zero base
    CIGAR = l[5]
    SEQ = l[9]
    #print QNAME, FLAG, RNAME, POS, CIGAR, SEQ
    rp = POS
    qp = 0
    for m in re.finditer(cigar_re, CIGAR):
        #print m.group(0)
        #print m.group(1), m.group(2)
        adv = int(m.group(1))
        #print m.group(0)
        if m.group(2) == "S":
            qp += adv
        if m.group(2) == "M":
            matches = []
            for i in range(adv):
                #print rp, qp, ref_seq[RNAME][rp], SEQ[qp]
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
                #print rp, qp, "-", SEQ[qp]
                qp += 1
        elif m.group(2) == "D":
            for i in range(adv):
                #print rp , qp, ref_seq[RNAME][rp], "-" 
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
                ref_base = ref_seq[RNAME][pos]
                print >> vpos2, pos+1, ref_base, total_count, " ".join(["%s %d" % (x[1], x[0]) for x in base_count])
                for q_id_ in pileup[pos][b0]:
                    print >> vmap2, pos+1, ref_base, b0, q_id_
                for q_id_ in pileup[pos][b1]:
                    print >> vmap2, pos+1, ref_base, b1, q_id_ 
            del pileup[pos]
    q_id += 1


q_id_map_f = open("q_id_map", "w")
for q_id, q_name in q_id_map.items():
    print >> q_id_map_f, q_id, q_name
