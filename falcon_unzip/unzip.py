from falcon_kit import run_support as support
#from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
#from pypeflow.task import PypeTask, PypeThreadTaskBase, PypeTaskBase
#from pypeflow.controller import PypeWorkflow, PypeThreadWorkflow
from pypeflow.simple_pwatcher_bridge import (
        PypeLocalFile, makePypeLocalFile, fn,
        PypeTask,
        PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase)
PypeThreadTaskBase = MyFakePypeThreadTaskBase
from falcon_kit.FastaReader import FastaReader
import glob
import logging
import os
import re
import sys
import time
import ConfigParser

LOG = logging.getLogger(__name__)

def system(call, check=False):
    LOG.debug('$(%s)' %repr(call))
    rc = os.system(call)
    msg = "Call %r returned %d." % (call, rc)
    if rc:
        LOG.warning(msg)
        if check:
            raise Exception(msg)
    else:
        LOG.debug(msg)
    return rc

def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def task_track_reads(self):
    job_done = fn(self.job_done)
    wd = self.parameters["wd"]
    config = self.parameters["config"]
    sge_track_reads = config["sge_track_reads"]
    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , "track_reads.sh")
    topdir = '../..'

    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date
cd {topdir}
python -m falcon_kit.mains.get_read_ctg_map
python -m falcon_kit.mains.rr_ctg_track
python -m falcon_kit.mains.pr_ctg_track
#mkdir -p 3-unzip/reads/
python -m falcon_kit.mains.fetch_reads
cd {wd}
date
touch {job_done}
""".format(**locals())

    with open(script_fn,"w") as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn
    #job_data["sge_option"] = sge_track_reads


def task_run_blasr(self):
    job_done = fn(self.job_done)
    ref_fasta = fn(self.ref_fasta)
    read_fasta = fn(self.read_fasta)

    job_uid = self.parameters["job_uid"]
    wd = self.parameters["wd"]
    ctg_id = self.parameters["ctg_id"]

    config = self.parameters["config"]
    smrt_bin = config["smrt_bin"]
    sge_blasr_aln = config["sge_blasr_aln"]
    job_type = config["job_type"]
    blasr = os.path.join(smrt_bin, "blasr")
    samtools = os.path.join( smrt_bin, "samtools")


    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , "aln_{ctg_id}.sh".format(ctg_id = ctg_id))

    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
cd {wd}
hostname
date
cd {wd}
time {blasr} {read_fasta} {ref_fasta} -noSplitSubreads -clipping subread\
 -hitPolicy randombest -randomSeed 42 -bestn 1 -minPctIdentity 70.0\
 -minMatch 12  -nproc 24 -sam -out tmp_aln.sam
{samtools} view -bS tmp_aln.sam | {samtools} sort - {ctg_id}_sorted
{samtools} index {ctg_id}_sorted.bam
rm tmp_aln.sam
date
touch {job_done}
""".format(**locals())

    with open(script_fn,"w") as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn
    #job_data["sge_option"] = sge_blasr_aln


def task_phasing(self):
    ref_fasta = fn(self.ref_fasta)
    aln_bam = fn(self.aln_bam)

    job_done = fn(self.job_done)

    job_uid = self.parameters["job_uid"]
    wd = self.parameters["wd"]
    ctg_id = self.parameters["ctg_id"]

    config = self.parameters["config"]
    sge_phasing = config["sge_phasing"]
    job_type = config["job_type"]

    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , "p_%s.sh" % (ctg_id))

    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
cd {wd}
hostname
date
cd {wd}
fc_phasing.py --bam {aln_bam} --fasta {ref_fasta} --ctg_id {ctg_id} --base_dir ../
fc_phasing_readmap.py --ctg_id {ctg_id} --read_map_dir ../../../2-asm-falcon/read_maps --phased_reads phased_reads
date
touch {job_done}
""".format(**locals())

    with open(script_fn,"w") as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn
    #job_data["sge_option"] = sge_phasing


def task_hasm(self):
    rid_to_phase_all = fn(self.rid_to_phase_all)
    job_done = fn(self.job_done)
    config = self.parameters["config"]
    sge_hasm = config["sge_hasm"]

    wd = self.parameters["wd"]

    job_type = config["job_type"]

    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , "hasm.sh")

    las_fofn = '../../2-asm-falcon/las.fofn'
    las_fofn = '../../1-preads_ovl/merge-gather/las.fofn'
    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date
cd {wd}

fc_ovlp_filter_with_phase.py --fofn {las_fofn} --max_diff 120 --max_cov 120 --min_cov 1 --n_core 12 --min_len 2500 --db ../../1-preads_ovl/preads.db --rid_phase_map {rid_to_phase_all} > preads.p_ovl
fc_phased_ovlp_to_graph.py preads.p_ovl --min_len 2500 > fc.log
if [ -e ../../1-preads_ovl/preads4falcon.fasta ];
then
  ln -sf ../../1-preads_ovl/preads4falcon.fasta .
else
  ln -sf ../../1-preads_ovl/db2falcon/preads4falcon.fasta .
fi
fc_graphs_to_h_tigs.py --fc_asm_path ../../2-asm-falcon/ --fc_hasm_path ./ --ctg_id all --rid_phase_map {rid_to_phase_all} --fasta preads4falcon.fasta

# more script -- a little bit hacky here, we should improve

WD=$PWD
for f in `cat ../reads/ctg_list `; do mkdir -p $WD/$f; cd $WD/$f; fc_dedup_h_tigs.py $f; done

## prepare for quviering the haplotig
cd $WD/..
if [ -e "all_phased_reads" ]; then rm all_phased_reads; fi
if [ -e "all_h_ctg_ids" ]; then rm all_h_ctg_ids; fi
if [ -e "all_p_ctg_edges" ]; then rm all_p_ctg_edges; fi
if [ -e "all_p_ctg.fa" ]; then rm all_p_ctg.fa; fi
if [ -e "all_h_ctg.fa" ]; then rm all_h_ctg.fa; fi

find 0-phasing -name "phased_reads" | sort | xargs cat >> all_phased_reads
find 1-hasm -name "h_ctg_ids.*" | sort | xargs cat >> all_h_ctg_ids
find 1-hasm -name "p_ctg_edges.*" | sort | xargs cat >> all_p_ctg_edges
find 1-hasm -name "h_ctg_edges.*" | sort | xargs cat >> all_h_ctg_edges
find 1-hasm -name "p_ctg.*.fa" | sort | xargs cat >> all_p_ctg.fa
find 1-hasm -name "h_ctg.*.fa" | sort | xargs cat >> all_h_ctg.fa
cd ../
date
touch {job_done}
""".format(**locals())

    with open(script_fn,"w") as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn
    #job_data["sge_option"] = sge_hasm

def unzip_all(config):
    unzip_concurrent_jobs = config["unzip_concurrent_jobs"]
    PypeProcWatcherWorkflow.setNumThreadAllowed(unzip_concurrent_jobs, unzip_concurrent_jobs)
    wf = PypeProcWatcherWorkflow()
    wf.max_jobs = unzip_concurrent_jobs

    ctg_list_file = makePypeLocalFile("./3-unzip/reads/ctg_list")
    falcon_asm_done = makePypeLocalFile("./2-asm-falcon/falcon_asm_done")
    wdir = os.path.abspath('./3-unzip/reads')
    parameters = {"wd": wdir, "config": config}

    job_done = makePypeLocalFile( os.path.join( parameters["wd"], "track_reads_done" ) )
    make_track_reads_task = PypeTask(inputs = {"falcon_asm_done": falcon_asm_done},
                                     outputs = {"job_done": job_done, "ctg_list_file": ctg_list_file},
                                     parameters = parameters,
                                     wdir = wdir,
                                     TaskType = PypeThreadTaskBase,
                                     URL = "task://localhost/track_reads" )
    track_reads_task = make_track_reads_task(task_track_reads)

    wf.addTask(track_reads_task)
    wf.refreshTargets() #force refresh now, will put proper dependence later

    ctg_ids = []
    with open("./3-unzip/reads/ctg_list") as f:
        for row in f:
            row = row.strip()
            ctg_ids.append( row )

    aln1_outs = {}

    all_ctg_out = {}

    for ctg_id in ctg_ids:
        # inputs
        ref_fasta = makePypeLocalFile("./3-unzip/reads/{ctg_id}_ref.fa".format(ctg_id = ctg_id))
        read_fasta = makePypeLocalFile("./3-unzip/reads/{ctg_id}_reads.fa".format(ctg_id = ctg_id))

        # outputs
        wd = os.path.join( os.getcwd(),  "./3-unzip/0-phasing/{ctg_id}/".format( ctg_id = ctg_id ) )
        #mkdir(wd)
        blasr_dir = os.path.join(wd, 'blasr')
        ctg_aln_out = makePypeLocalFile( os.path.join( blasr_dir, "{ctg_id}_sorted.bam".format( ctg_id = ctg_id ) ) )
        job_done = makePypeLocalFile( os.path.join( blasr_dir, "aln_{ctg_id}_done".format( ctg_id = ctg_id ) ) )

        parameters = {"job_uid":"aln-"+ctg_id, "wd": blasr_dir, "config":config, "ctg_id": ctg_id}
        make_blasr_task = PypeTask(inputs = {"ref_fasta": ref_fasta, "read_fasta": read_fasta},
                                   outputs = {"ctg_aln_out": ctg_aln_out, "job_done": job_done},
                                   parameters = parameters,
                                   TaskType = PypeThreadTaskBase,
                                   URL = "task://localhost/aln_{ctg_id}".format( ctg_id = ctg_id ) )
        blasr_task = make_blasr_task(task_run_blasr)
        aln1_outs[ctg_id] = (ctg_aln_out, job_done)
        wf.addTask(blasr_task)

        phasing_dir = os.path.join(wd, 'phasing')
        job_done = makePypeLocalFile( os.path.join( phasing_dir, "p_{ctg_id}_done".format( ctg_id = ctg_id ) ) )
        rid_to_phase_out = makePypeLocalFile( os.path.join( wd, "rid_to_phase.{ctg_id}".format( ctg_id = ctg_id ) ) ) # TODO: ???
        all_ctg_out[ "r2p.{ctg_id}".format( ctg_id = ctg_id ) ] = rid_to_phase_out # implicit output?

        parameters = {"job_uid":"ha-"+ctg_id, "wd": wd, "config":config, "ctg_id": ctg_id}
        make_phasing_task = PypeTask(inputs = {"ref_fasta": ref_fasta, "aln_bam":ctg_aln_out},
                                   outputs = {"job_done": job_done},
                                   parameters = parameters,
                                   TaskType = PypeThreadTaskBase,
                                   URL = "task://localhost/p_{ctg_id}".format( ctg_id = ctg_id ) )
        phasing_task = make_phasing_task(task_phasing)
        wf.addTask(phasing_task)

    wf.refreshTargets()

    hasm_wd = os.path.abspath("./3-unzip/1-hasm/")
    #mkdir(hasm_wd)
    rid_to_phase_all =  makePypeLocalFile( os.path.join(hasm_wd, 'rid-to-phase-all', "rid_to_phase.all") )
    task = PypeTask(inputs = all_ctg_out, outputs = {"rid_to_phase_all": rid_to_phase_all},
                TaskType = PypeThreadTaskBase, URL = "task://localhost/rid_to_phase_all" )(get_rid_to_phase_all)
    wf.addTask(task)

    parameters["wd"] = hasm_wd
    job_done = makePypeLocalFile( os.path.join( hasm_wd, "hasm_done" ) )
    make_hasm_task = PypeTask(inputs = {"rid_to_phase_all": rid_to_phase_all},
                              outputs = {"job_done": job_done},
                              parameters = parameters,
                              TaskType = PypeThreadTaskBase,
                              URL = "task://localhost/hasm" )
    hasm_task = make_hasm_task(task_hasm)

    wf.addTask(hasm_task)

    wf.refreshTargets()

def get_rid_to_phase_all(self):
    # Tasks must be at module scope now.
    rid_to_phase_all_fn = fn(self.rid_to_phase_all)
    inputs_fn = [ fn(f) for f in self.inputs.values() ]
    inputs_fn.sort()
    output = []
    for fname in inputs_fn:
        output.extend( open(fname).read() )

    out = open( rid_to_phase_all_fn, "w")
    out.write("".join(output))
    out.close()

def main(argv=sys.argv):

    if len(argv) < 2 or argv[1].startswith('-'):
        print "you need to provide a configuration file to specific a couple cluster running environment"
        sys.exit(1)

    config_fn = argv[1]

    config = ConfigParser.ConfigParser()
    config.read(config_fn)

    job_type = "SGE"
    if config.has_option('General', 'job_type'):
        job_type = config.get('General', 'job_type')

    sge_blasr_aln = " -pe smp 24 -q bigmem "
    if config.has_option('Unzip', 'sge_blasr_aln'):
        sge_blasr_aln = config.get('Unzip', 'sge_blasr_aln')

    smrt_bin = "/mnt/secondary/builds/full/3.0.0/prod/smrtanalysis_3.0.0.153854/smrtcmds/bin/"
    if config.has_option('Unzip', 'smrt_bin'):
        smrt_bin = config.get('Unzip', 'smrt_bin')

    sge_phasing = " -pe smp 12 -q bigmem"
    if config.has_option('Unzip', 'sge_phasing'):
        sge_phasing = config.get('Unzip', 'sge_phasing')

    sge_hasm = " -pe smp 48 -q bigmem"
    if config.has_option('Unzip', 'sge_hasm'):
        sge_hasm = config.get('Unzip', 'sge_hasm')

    sge_track_reads = " -pe smp 12 -q bigmem"
    if config.has_option('Unzip', 'sge_track_reads'):
        sge_track_reads = config.get('Unzip', 'sge_track_reads')

    unzip_concurrent_jobs = 8
    if config.has_option('Unzip', 'unzip_concurrent_jobs'):
        unzip_concurrent_jobs = config.getint('Unzip', 'unzip_concurrent_jobs')

    config = {"job_type": job_type,
              "sge_blasr_aln": sge_blasr_aln,
              "smrt_bin": smrt_bin,
              "sge_phasing": sge_phasing,
              "sge_hasm": sge_hasm,
              "sge_track_reads": sge_track_reads,
              "unzip_concurrent_jobs":unzip_concurrent_jobs}

    #support.job_type = "SGE" #tmp hack until we have a configuration parser

    unzip_all(config)
