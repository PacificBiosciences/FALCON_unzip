from falcon_kit import run_support as support
from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
from pypeflow.task import PypeTask, PypeThreadTaskBase, PypeTaskBase
from pypeflow.controller import PypeWorkflow, PypeThreadWorkflow
from falcon_kit.FastaReader import FastaReader
import glob
import os
import re
import sys
import time
import ConfigParser


global fc_run_logger
fc_run_logger = support.setup_logger(None)

support.job_type = "SGE" #tmp hack until we have a configuration parser

wait_time = 5
fc_run_logger = None

def system(call, check=False):
    fc_run_logger.debug('$(%s)' %repr(call))
    rc = os.system(call)
    msg = "Call %r returned %d." % (call, rc)
    if rc:
        fc_run_logger.warning(msg)
        if check:
            raise Exception(msg)
    else:
        fc_run_logger.debug(msg)
    return rc

def _qsub_script(job_data, specific):
        script_fn = job_data["script_fn"]
        job_name = job_data["job_name"]
        cwd = job_data["cwd"]
        sge_option = job_data["sge_option"]
        sge_cmd="qsub -N {job_name} {sge_option} -o {cwd}/sge_log {specific}\
                 -S /bin/bash {script}".format(job_name=job_name,
                                               cwd=os.getcwd(),
                                               specific=specific,
                                               sge_option=sge_option,
                                               script=script_fn)
        system(sge_cmd, check=True)
        system("sleep 1")

def _run_script_sge(job_data):
    specific = '-j y'
    _qsub_script(job_data, specific)

def _run_script_torque(job_data):
    # See https://github.com/PacificBiosciences/FALCON/pull/227
    specific = '-j oe'
    _qsub_script(job_data, specific)

def _run_script_slurm(job_data):
        script_fn = job_data["script_fn"]
        job_name = job_data["job_name"]
        cwd = job_data["cwd"]
        sge_option = job_data["sge_option"]
        with open(script_fn, 'r') as original: data = original.read()
        with open(script_fn, 'w') as modified: modified.write("#!/bin/sh" + "\n" + data)
        sge_cmd="sbatch -J {job_name} {sge_option} {script}".format(job_name=job_name, cwd=os.getcwd(),sge_option=sge_option, script=script_fn)
        system(sge_cmd, check=True)

def _run_script_local(job_data):
        script_fn = job_data["script_fn"]
        job_name = job_data["job_name"]
        log_fn = '{0}.log'.format(script_fn)
        cmd = "bash {0} 1> {1} 2>&1".format(script_fn, log_fn)
        try:
            system(cmd, check=True)
        except Exception:
            out = open(log_fn).read()
            fc_run_logger.exception('Contents of %r:\n%s' %(log_fn, out))
            raise

_run_scripts = {
        'SGE': _run_script_sge,
        'TORQUE': _run_script_torque,
        'SLURM': _run_script_slurm,
        'LOCAL': _run_script_local,
}

def run_script(job_data, job_type = "SGE" ):
    """For now, we actually modify the script before running it.
    This assume a simple bash script.
    We will have a better solution eventually.
    """
    try:
        _run_script = _run_scripts[job_type.upper()]
    except LookupError as e:
        msg = 'Unknown job_type=%s' %repr(job_type)
        fc_run_logger.exception(msg)
        raise
    job_name = job_data["job_name"]
    script_fn = job_data["script_fn"]
    support.update_env_in_script(script_fn,
        ['PATH', 'PYTHONPATH', 'LD_LIBRARY_PATH'])
    fc_run_logger.info('(%s) %r' %(job_type, script_fn))
    fc_run_logger.debug('%s (job %r)' %(_run_script.__name__, job_name))
    rc = _run_script(job_data)
    # Someday, we might trap exceptions here, as a failure would be caught later anyway.

def mkdir(d):
    if not os.path.isdir(d):
        os.makedirs(d)

def wait_for_file(filename, task, job_name = ""):
    """We could be in the thread or sub-process which spawned a qsub job,
    so we must check for the shutdown_event.
    """
    while 1:
        time.sleep(wait_time)
        # We prefer all jobs to rely on `*done.exit`, but not all do yet. So we check that 1st.
        exit_fn = filename + '.exit'
        if os.path.exists(exit_fn):
            fc_run_logger.info( "%r found." % (exit_fn) )
            fc_run_logger.debug( " job: %r exited." % (job_name) )
            os.unlink(exit_fn) # to allow a restart later, if not done
            if not os.path.exists(filename):
                fc_run_logger.warning( "%r is missing. job: %r failed!" % (filename, job_name) )
            break
        if os.path.exists(filename) and not os.path.exists(exit_fn):
            # (rechecked exit_fn to avoid race condition)
            fc_run_logger.info( "%r not found, but job is done." % (exit_fn) )
            fc_run_logger.debug( " job: %r exited." % (job_name) )
            break
        if task.shutdown_event is not None and task.shutdown_event.is_set():
            fc_run_logger.warning( "shutdown_event received (Keyboard Interrupt maybe?), %r not finished."
                % (job_name) )
            if support.job_type == "SGE":
                fc_run_logger.info( "deleting the job by `qdel` now..." )
                system("qdel %s" % job_name) # Failure is ok.
            if support.job_type == "SLURM":
                fc_run_logger.info( "Deleting the job by 'scancel' now...")
                system("scancel -n %s" % job_name)
            break

def task_track_reads(self):

    job_done = fn(self.job_done)
    wd = self.parameters["wd"]
    config = self.parameters["config"]
    input_bam_fofn = config["input_bam_fofn"]
    sge_track_reads = config["sge_track_reads"]
    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , "track_reads_h.sh")

    script = []
    script.append( "set -vex" )
    script.append( "trap 'touch {job_done}.exit' EXIT".format(job_done = job_done) )
    script.append( "cd %s" % wd )
    script.append( "hostname" )
    script.append( "date" )
    script.append( "cd {wd}".format(wd = wd) )
    script.append( "fc_get_read_hctg_map.py" )
    script.append( "fc_rr_hctg_track.py" )
    script.append( "mkdir -p 4-quiver/reads/" )
    script.append( "fc_select_reads_from_bam.py {input_bam_fofn}".format( input_bam_fofn = input_bam_fofn) )
    script.append( "date" )
    script.append( "touch {job_done}".format(job_done = job_done) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script) + '\n')

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_track_reads
    run_script(job_data, job_type = config["job_type"])
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])


def task_run_quiver(self):

    ref_fasta = fn(self.ref_fasta)
    read_sam = fn(self.read_sam)

    cns_fasta = fn(self.cns_fasta)
    cns_fastq = fn(self.cns_fastq)
    job_done = fn(self.job_done)

    job_uid = self.parameters["job_uid"]
    wd = self.parameters["wd"]
    config = self.parameters["config"]
    ctg_id = self.parameters["ctg_id"]

    smrt_bin = config["smrt_bin"]
    sge_quiver = config["sge_quiver"]
    job_type = config["job_type"]
    samtools = os.path.join( smrt_bin, "samtools")
    pbalign = os.path.join( smrt_bin, "pbalign")
    makePbi = os.path.join( smrt_bin, "makePbi")
    variantCaller = os.path.join( smrt_bin, "variantCaller")

    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , "cns_%s.sh" % (ctg_id))

    script = []
    script.append( "set -vex" )
    script.append( "trap 'touch {job_done}.exit' EXIT".format(job_done = job_done) )
    script.append( "cd %s" % wd )
    script.append( "hostname" )
    script.append( "date" )
    script.append( "cd {wd}".format(wd = wd) )

    script.append( "{samtools} faidx {ref_fasta}".format( samtools=samtools, ref_fasta=ref_fasta ) )
    script.append( "{samtools} view -b -S {read_sam} > {ctg_id}.bam".format( samtools=samtools, read_sam = read_sam, ctg_id = ctg_id ) )
    script.append( "{pbalign} --tmpDir=/localdisk/scratch/ --nproc=24 --minAccuracy=0.75 --minLength=50\
            --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr\
            --algorithmOptions=-useQuality --maxHits=1 --hitPolicy=random --seed=1\
            {ctg_id}.bam {ref_fasta} aln-{ctg_id}.bam".format( pbalign=pbalign , ctg_id = ctg_id, ref_fasta = ref_fasta))
    script.append( "#{makePbi} --referenceFasta {ref_fasta} aln-{ctg_id}.bam".format(makePbi = makePbi, ref_fasta = ref_fasta, ctg_id = ctg_id) )
    script.append( "({variantCaller} -x 5 -X 120 -q 20 -j 24 -r {ref_fasta} aln-{ctg_id}.bam\
            -o {cns_fasta} -o {cns_fastq}) || echo quvier failed".format( variantCaller = variantCaller, ctg_id = ctg_id, ref_fasta = ref_fasta,
                                                   cns_fasta=cns_fasta, cns_fastq=cns_fastq ))

    script.append( "date" )
    script.append( "touch {job_done}".format(job_done = job_done) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script) + '\n')

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_quiver
    run_script(job_data, job_type = job_type)
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])


def main(argv=sys.argv):

    global fc_run_logger
    fc_run_logger = support.setup_logger(None)


    if len(sys.argv) < 2:
        print "you need to provide a configuration file to specific a couple cluster running environment"
        sys.exit(1)

    config_fn = sys.argv[1]

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


    quiver_concurrent_jobs = 8
    if config.has_option('Unzip', 'quiver_concurrent_jobs'):
        quiver_concurrent_jobs = config.getint('Unzip', 'quiver_concurrent_jobs')

    config = {"job_type": job_type,
              "sge_quiver": sge_quiver,
              "sge_track_reads": sge_track_reads,
              "input_bam_fofn": input_bam_fofn,
              "smrt_bin": smrt_bin}

    support.job_type = "SGE" #tmp hack until we have a configuration parser

    ctg_ids = []


    PypeThreadWorkflow.setNumThreadAllowed(quiver_concurrent_jobs, quiver_concurrent_jobs)
    wf = PypeThreadWorkflow()

    parameters = {"wd": os.path.abspath("."), "config": config}
    hasm_done = makePypeLocalFile("./3-unzip/1-hasm/hasm_done")
    job_done = makePypeLocalFile( os.path.join( parameters["wd"], "track_reads_h_done" ) )
    make_track_reads_task = PypeTask(inputs = {"hasm_done": hasm_done},
                                     outputs = {"job_done": job_done},
                                     parameters = parameters,
                                     TaskType = PypeThreadTaskBase,
                                     URL = "task://localhost/track_reads_h" )
    track_reads_task = make_track_reads_task(task_track_reads)

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



    wf.refreshTargets()
    os.system("sleep 30")

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
