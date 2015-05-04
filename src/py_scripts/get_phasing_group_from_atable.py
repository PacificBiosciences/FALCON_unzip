

left_connect = {}
right_connect = {}

c_score = {}
states = {}
positions = set() 

def get_score( c_score, pos1, pos2, s1, s2 ):
    if pos1 > pos2:
        pos1, pos2 = pos2, pos1
        s1, s2 = s2, s1
    b11, b12 = s1
    b21, b22 = s2
    return c_score[ (pos1, pos2) ][ (b11+b21, b12+b22) ]


ref_base = {}
with open("vmap2") as f:
    for l in f:
        l = l.strip().split()
        pos = int(l[0])
        ref_b = l[1]
        v_b = l[2]
        q_id = int(l[3])
        ref_base[pos] = ref_b

with open("atable") as f:
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
            connects = left_connect.get(pos1,[]) + right_connect.get(pos1,[]) 
            if len(connects) == 0:
                states[pos1] = (b11, b12)
            else:
                st1 = (b11, b12)
                st2 = (b12, b11)
                score1 = 0
                score2 = 0
                for pp in connects:
                    if pp in states:
                        st0 = states[pp]
                    else:
                        continue
                    score1 += get_score( c_score, pos1, pp, st1, st0 )
                    score2 += get_score( c_score, pos1, pp, st2, st0 )
                if score1 > score2:
                    states[pos1] = st1
                else:
                    states[pos1] = st2


        if pos2 not in states:
            connects = left_connect.get(pos2,[]) + right_connect.get(pos2,[]) 
            if len(connects) == 0:
                states[pos2] = (b21, b22)
            else:
                st1 = (b21, b22)
                st2 = (b22, b21)
                score1 = 0
                score2 = 0
                for pp in connects:
                    if pp in states:
                        st0 = states[pp]
                    else:
                        continue
                    score1 += get_score( c_score, pp, pos2, st0, st1 )
                    score2 += get_score( c_score, pp, pos2, st0, st2 )
                if score1 > score2:
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
        for pp in left_connect.get(p,[]) + right_connect.get(p,[]):
            st0 = states[pp]
            score1 += get_score( c_score, pp, p, st0 ,st1)
            score2 += get_score( c_score, pp, p, st0, st2)

        if score1 >= score2:
            states[p] = st1
        else:
            states[p] = st2
            update_count += 1
    if update_count == 0:
        break


right_extent = {}
left_extent = {}

for p in positions:

    left_extent[p] = p
    if p in left_connect:
        left = p
        st0 = states[p]
        for pp in left_connect[p]:
            st1 = states[pp]
            if get_score( c_score, p, pp, st0, st1) > 0 and pp < left:
                left = pp
        left_extent[p] = left
            
    right_extent[p] = p
    if p in right_connect:
        right = p
        st0 = states[p]
        for pp in right_connect[p]:
            st1 = states[pp]
            if get_score( c_score, p, pp, st0, st1) > 0 and pp > right:
                right = pp
        right_extent[p] = right


phase_block_id = 1
phase_blocks = {}
pb = []
"""
adv = 0
for p in positions:
    b1, b2 = states[p]

    if len(right_connect.get(p,[])) > 0 or adv > p:
        pb.append( (p, b1, b2) )
        if right_extent[p] > adv:
            adv = right_extent[p] 
    else:
        if len(pb) > 1:
            phase_blocks[phase_block_id] = pb
            phase_block_id += 1
        pb = []
    if len(left_connect.get(p,[])) == 0 and adv < p:
        pb = []

phase_blocks[phase_block_id] = pb
"""

max_right_ext = 0
for p in positions:
    b1, b2 = states[p]
    if max_right_ext < left_extent[p]:
        if len(pb) > 0:
            phase_blocks[phase_block_id] = pb
            phase_block_id += 1
        pb = []
    pb.append( (p, b1, b2) )
    if right_extent[p] > max_right_ext:
        max_right_ext =  right_extent[p]

phase_blocks[phase_block_id] = pb

for pid in xrange(1, phase_block_id+1):
    if len(phase_blocks[pid]) == 0:
        continue
    min_ = min( [x[0] for x in phase_blocks[pid]] )
    max_ = max( [x[0] for x in phase_blocks[pid]] )
    
    print "P", pid, min_, max_, max_ - min_, len(phase_blocks[pid]), 1.0 * (max_-min_)/len(phase_blocks[pid]) 
    for p, b1, b2 in phase_blocks[pid]:
        rb = ref_base[p]
        print "V", pid, p, "%d_%s_%s" % (p,rb,b1), "%d_%s_%s" % (p,rb,b2), left_extent[p], right_extent[p] 
