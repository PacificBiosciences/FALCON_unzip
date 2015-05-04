import sys

mum_fn = sys.argv[1]

def process_mum_data( pread_id, mum_data ):
    delta = 0
    pend = 0
    qend = 0
    p_delta = 0


    rtn = []
    for rpos, qpos, len_ in mum_data:
        rpos = int(rpos)
        qpos = int(qpos)
        len_ = int(len_)
        offset = rpos - qpos
        
        if pend == 0:
            rtn.append( ( rpos, qpos, len_, offset, rpos+len_, 0, 0, 0, pread_id) )
        else:
            rtn.append( ( rpos, qpos, len_, offset, rpos+len_, offset-p_offset, rpos-pend, qpos-qend, pread_id) )
        pend = rpos+len_
        qend = qpos+len_
        p_offset = offset

    #break into cluster
    cluster = []
    p_delta = 0 
    for d in rtn:
        rpos, qpos, len_, offset, r_end, delta, r_delta, q_delta, pread_id = d
        if q_delta < -10:
            continue
        if r_delta < -10:
            continue
        if len(cluster) == 0: #start a new cluster
            cluster.append([])
            cluster[-1].append( d )
            #p_delta = delta
            continue
        if delta  > 30000:
            #print delta, p_delta
            #p_delta = 0
            cluster.append([])
            continue
        cluster[-1].append( d )
        #p_delta = delta
    rtn = []
    for cls_id in xrange(len(cluster)):
        aln_size = 0
        qpos_ary = []
        if len(cluster[cls_id]) == 0:
            continue
        for d in cluster[cls_id]:
            rpos, qpos, len_, offset, r_end, delta, r_delta, q_delta, pread_id = d
            aln_size += len_
            qpos_ary.append(qpos)
        if 1.0*(max(qpos_ary) - min(qpos_ary))/aln_size < 0.5 or aln_size < 2000:
            continue

        for d in cluster[cls_id]:
            rpos, qpos, len_, offset, r_end, delta, r_delta, q_delta, pread_id = d
            rtn.append( (rpos, qpos, len_, offset, r_end, delta, r_delta, q_delta, pread_id) )


    return rtn



mum_data = []
with open(mum_fn) as mum_file:
    for row in mum_file:
        row = row.strip().split()
        if row[0] == ">":
            if len(mum_data) > 0:
                rtn = process_mum_data( pread_id, mum_data )
                if len(rtn) > 0:
                    print ">",pread_id
                    for d in rtn:
                        print " ".join([str(c) for c in d])
            pread_id = row[1]
            mum_data = []
            #print 
        else:
            mum_data.append( row )
