import subprocess, shlex


read_to_phase = {}
with open("phased_reads") as f:
    for l in f:
        l = l.strip().split()
        rid = l[-1]
        phase = l[1] + "_" +l[2]
        read_to_phase[rid] = phase

p = subprocess.Popen(shlex.split("samtools view -h 000207F.000207F000207F.bam" ), stdout=subprocess.PIPE)

phase_to_file = {}
header = []
for l in p.stdout:
    l2 = l.strip().split()
    if l2[0][0] == "@":
        header.append(l.strip())
    QNAME = l2[0]
    phase = read_to_phase.get(QNAME, None)
    if phase != None:
        pfn = ("./phased_sam/%s" % phase) +".sam"
        if pfn not in phase_to_file:
            phase_to_file[pfn] = open(pfn,"w")
            print >> phase_to_file[pfn], "\n".join(header)
        print >>  phase_to_file[pfn], l.strip()

