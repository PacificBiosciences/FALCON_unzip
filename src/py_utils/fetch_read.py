from falcon_kit.FastaReader import FastaReader

import sys

f = FastaReader(sys.argv[1])
rl = set(open(sys.argv[2]).read().split())
for r in f:
    rid = r.name.split()[0]
    if rid not in rl:
        continue
    print ">"+rid
    print r.sequence
