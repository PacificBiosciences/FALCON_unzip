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
import pprint
import re
import sys
import time
import ConfigParser


LOG = logging.getLogger(__name__)

def system(call, check=False):
    LOG.debug('$(%s)' %repr(call))
    rc = os.system(call)
    msg = 'Call %r returned %d.' % (call, rc)
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
    wd = self.parameters['wd']
    config = self.parameters['config']
    input_bam_fofn = config['input_bam_fofn']
    sge_track_reads = config['sge_track_reads']
    script_dir = os.path.join(wd)
    script_fn = os.path.join(script_dir, 'track_reads_h.sh')

    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
cd {wd}
hostname
date
fc_get_read_hctg_map.py --basedir ../..
fc_rr_hctg_track.py --base_dir ../..
mkdir -p 4-quiver/reads/
fc_select_reads_from_bam.py --basedir ../.. {input_bam_fofn}
date
touch {job_done}
""".format(**locals())

    with open(script_fn,'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn
    #job_data["sge_option"] = sge_track_reads


def task_run_quiver(self):

    ref_fasta = fn(self.ref_fasta)
    read_sam = fn(self.read_sam)

    cns_fasta = fn(self.cns_fasta)
    cns_fastq = fn(self.cns_fastq)
    job_done = fn(self.job_done)

    job_uid = self.parameters['job_uid']
    wd = self.parameters['wd']
    config = self.parameters['config']
    ctg_id = self.parameters['ctg_id']

    smrt_bin = config['smrt_bin']
    sge_quiver = config['sge_quiver']
    job_type = config['job_type']
    samtools = os.path.join( smrt_bin, 'samtools')
    pbalign = os.path.join( smrt_bin, 'pbalign')
    makePbi = os.path.join( smrt_bin, 'makePbi')
    variantCaller = os.path.join( smrt_bin, 'variantCaller')

    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , 'cns_%s.sh' % (ctg_id))

    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date
cd {wd}

{samtools} faidx {ref_fasta}
{samtools} view -b -S {read_sam} > {ctg_id}.bam
{pbalign} --tmpDir=/localdisk/scratch/ --nproc=24 --minAccuracy=0.75 --minLength=50\
          --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr\
          --algorithmOptions=-useQuality --maxHits=1 --hitPolicy=random --seed=1\
            {ctg_id}.bam {ref_fasta} aln-{ctg_id}.bam
#{makePbi} --referenceFasta {ref_fasta} aln-{ctg_id}.bam
({variantCaller} -x 5 -X 120 -q 20 -j 24 -r {ref_fasta} aln-{ctg_id}.bam\
            -o {cns_fasta} -o {cns_fastq}) || echo quvier failed
date
touch {job_done}
""".format(**locals())

    with open(script_fn,'w') as script_file:
        script_file.write(script)
    self.generated_script_fn = script_fn
    #job_data["sge_option"] = sge_quiver


def main(argv=sys.argv):
    global LOG
    LOG = support.setup_logger(None)


    if len(sys.argv) < 2:
        print "you need to provide a configuration file to specific a couple cluster running environment"
        sys.exit(1)

    config_fn = sys.argv[1]
    config_absbasedir = os.path.dirname(os.path.abspath(config_fn))

    config = ConfigParser.ConfigParser()
    config.read(config_fn)


    job_type = "SGE"
    if config.has_option('General', 'job_type'):
        job_type = config.get('General', 'job_type')

    sge_track_reads = " -pe smp 12 -q bigmem"
    if config.has_option('Unzip', 'sge_track_reads'):
        sge_track_reads = config.get('Unzip', 'sge_track_reads')

    sge_quiver = " -pe smp 24 -q bigmem "
    if config.has_option('Unzip', 'sge_quiver'):
        sge_quiver = config.get('Unzip', 'sge_quiver')

    smrt_bin = "/mnt/secondary/builds/full/3.0.0/prod/smrtanalysis_3.0.0.153854/smrtcmds/bin/"
    if config.has_option('Unzip', 'smrt_bin'):
        smrt_bin = config.get('Unzip', 'smrt_bin')

    input_bam_fofn = "input_bam.fofn"
    if config.has_option('Unzip', 'input_bam_fofn'):
        input_bam_fofn = config.get('Unzip', 'input_bam_fofn')
    if not os.path.isabs(input_bam_fofn):
        input_bam_fofn = os.path.join(config_absbasedir, input_bam_fofn)


    quiver_concurrent_jobs = 8
    if config.has_option('Unzip', 'quiver_concurrent_jobs'):
        quiver_concurrent_jobs = config.getint('Unzip', 'quiver_concurrent_jobs')

    config = {"job_type": job_type,
              "sge_quiver": sge_quiver,
              "sge_track_reads": sge_track_reads,
              "input_bam_fofn": input_bam_fofn,
              "smrt_bin": smrt_bin}
    LOG.info('config={}'.format(pprint.pformat(config)))

    support.job_type = "SGE" #tmp hack until we have a configuration parser

    ctg_ids = []


    PypeProcWatcherWorkflow.setNumThreadAllowed(quiver_concurrent_jobs, quiver_concurrent_jobs)
    wf = PypeProcWatcherWorkflow()
    wf.max_jobs = quiver_concurrent_jobs

    abscwd = os.path.abspath('.')
    parameters = {'wd': os.path.join(abscwd, '4-quiver', 'track_reads_h'), 'config': config}
    hasm_done = makePypeLocalFile('./3-unzip/1-hasm/hasm_done')
    job_done = makePypeLocalFile( os.path.join( parameters['wd'], 'track_reads_h_done' ) )
    make_track_reads_task = PypeTask(inputs = {'hasm_done': hasm_done},
                                     outputs = {'job_done': job_done},
                                     parameters = parameters,
                                     TaskType = PypeThreadTaskBase,
                                     URL = 'task://localhost/track_reads_h' )
    track_reads_task = make_track_reads_task(task_track_reads)
    #sge_track_reads = config["sge_track_reads"]

    wf.addTask(track_reads_task)
    wf.refreshTargets() #force refresh now, will put proper dependence later

    ref_seq_data = {}
    p_ctg_fa = FastaReader("./3-unzip/all_p_ctg.fa")
    ctg_types = {}
    for r in p_ctg_fa:
        rid = r.name.split()[0]
        ref_seq_data[rid] = r.sequence
        ctg_types[rid] = "p"


    h_ctg_fa = FastaReader("./3-unzip/all_h_ctg.fa")
    for r in h_ctg_fa:
        rid = r.name.split()[0]
        ref_seq_data[rid] = r.sequence
        ctg_types[rid] = "h"

    ctg_ids = sorted( ref_seq_data.keys() )
    p_ctg_out=[]
    h_ctg_out=[]
    for ctg_id in ctg_ids:
        sequence = ref_seq_data[ctg_id]
        m_ctg_id = ctg_id.split("-")[0]
        wd = os.path.join( os.getcwd(), "./4-quiver/", m_ctg_id )
        mkdir( wd )
        ref_fasta = makePypeLocalFile(os.path.join(wd, "{ctg_id}_ref.fa".format(ctg_id = ctg_id) ) )
        read_sam = makePypeLocalFile(os.path.join( os.getcwd(), "./4-quiver/reads/" "{ctg_id}.sam".format(ctg_id = ctg_id) ) )
        cns_fasta = makePypeLocalFile(os.path.join(wd, "cns-{ctg_id}.fasta.gz".format(ctg_id = ctg_id) ) )
        cns_fastq = makePypeLocalFile(os.path.join(wd, "cns-{ctg_id}.fastq.gz".format(ctg_id = ctg_id) ) )
        job_done = makePypeLocalFile(os.path.join(wd, "{ctg_id}_quiver_done".format(ctg_id = ctg_id) ) )

        if os.path.exists(fn(read_sam)):
            if ctg_types[ctg_id] == "p":
                p_ctg_out.append( (cns_fasta, cns_fastq) )
            if ctg_types[ctg_id] == "h":
                h_ctg_out.append( (cns_fasta, cns_fastq) )
            if not os.path.exists(fn(ref_fasta)):
                with open(fn(ref_fasta),"w") as f:
                    print >>f, ">"+ctg_id
                    print >>f, sequence
            parameters = {"job_uid":"q-"+ctg_id, "wd": wd, "config":config, "ctg_id": ctg_id }
            make_quiver_task = PypeTask(inputs = {"ref_fasta": ref_fasta, "read_sam": read_sam},
                                       outputs = {"cns_fasta": cns_fasta, "cns_fastq": cns_fastq, "job_done": job_done},
                                       parameters = parameters,
                                       TaskType = PypeThreadTaskBase,
                                       URL = "task://localhost/q_{ctg_id}".format( ctg_id = ctg_id ) )
            quiver_task = make_quiver_task(task_run_quiver)
            wf.addTask( quiver_task )
    #sge_quiver = config["sge_quiver"]


    wf.refreshTargets()
    #os.system("sleep 30")

    mkdir( "./4-quiver/cns_output" )
    os.system("rm ./4-quiver/cns_output/cns_p_ctg.fasta")
    os.system("rm ./4-quiver/cns_output/cns_p_ctg.fastq")
    for cns_fasta, cns_fastq in sorted(p_ctg_out):
        os.system("zcat {cns_fasta} >> ./4-quiver/cns_output/cns_p_ctg.fasta".format( cns_fasta=fn(cns_fasta) ) )
        os.system("zcat {cns_fastq} >> ./4-quiver/cns_output/cns_p_ctg.fastq".format( cns_fastq=fn(cns_fastq) ) )



    os.system("rm ./4-quiver/cns_output/cns_h_ctg.fasta")
    os.system("rm ./4-quiver/cns_output/cns_h_ctg.fastq")
    for cns_fasta, cns_fastq in sorted(h_ctg_out):
        os.system("zcat {cns_fasta} >> ./4-quiver/cns_output/cns_h_ctg.fasta".format( cns_fasta=fn(cns_fasta) ) )
        os.system("zcat {cns_fastq} >> ./4-quiver/cns_output/cns_h_ctg.fastq".format( cns_fastq=fn(cns_fastq) ) )
