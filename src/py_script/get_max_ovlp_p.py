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
    for row in sp.check_output(shlex.split("LA4Falcon -m preads.db preads.%d.las" % c_id) ).splitlines():
        row = row.strip().split()
        print >>f, " ".join(row)
