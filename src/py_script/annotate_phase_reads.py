from collections import Counter
mn_to_sample = {}

with open("CHM13_norm_ids") as f:
    for row in f:
        row = row.strip().split()
        mn = row[1].split("/")[0]
        if mn not in mn_to_sample:
            mn_to_sample[mn] = "CHM13"

with open("CHM1_norm_ids") as f:
    for row in f:
        row = row.strip().split()
        mn = row[1].split("/")[0]
        if mn not in mn_to_sample:
            mn_to_sample[mn] = "CHM1"

phase_count = {}
with open("phased_reads") as f:
    for row in f:
        row = row.strip()
        row_s = row.split()
        rid = row_s[5]
        mn = rid.split("/")[0]
        phasing_block = int(row_s[1]) 
        phase = int(row_s[2]) 
        phase_count.setdefault( phasing_block, {0:[], 1:[]} )

        phase_count[phasing_block][phase].append( mn_to_sample[mn] )
        print "R", row, mn_to_sample[mn] 

for phasing_block in phase_count:
    cnt = Counter(phase_count[phasing_block][0] )
    print "P", phasing_block, 0, cnt["CHM1"], cnt["CHM13"]
    cnt = Counter(phase_count[phasing_block][1] )
    print "P", phasing_block, 1, cnt["CHM1"], cnt["CHM13"]
