import sys
import subprocess as sp
import shlex

c_id = int(sys.argv[1])
a_id = None
b_id = None
longest_ovlp_data = []
ovlp_count = 0
longest_ovlp = 0
with open("max_ovlp.%d" % c_id, "w") as f:
    for row in sp.check_output(shlex.split("LA4Falcon -mo raw_reads.db raw_reads.%d.las " % c_id) ).splitlines():
        row = row.strip().split()
        if row[-1] == "contained":
            continue

        if row[0] != a_id:
            if a_id != None and len(longest_ovlp_data) > 0:
                longest_ovlp_data.sort()
                for d in longest_ovlp_data[:200]:
                    print >>f, " ".join(d[1]), ovlp_count
            ovlp_count = 0
            longest_ovlp = 0
            longest_ovlp_data = []
            a_id = row[0]

        ovlp_len  = int(row[6]) - int(row[5])
        longest_ovlp_data.append( (-ovlp_len, row) )
        ovlp_count += 1


