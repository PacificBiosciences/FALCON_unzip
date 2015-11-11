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



def task_run_blasr(self):

    job_done = fn(self.job_done)
    ref_fasta = fn(self.ref_fasta)
    read_fasta = fn(self.read_fasta)

    job_uid = self.parameters["job_uid"]
    wd = self.parameters["wd"]
    sge_blasr_aln = self.parameters["sge_blasr_aln"]
    ctg_id = self.parameters["ctg_id"]
    smrt_bin = self.parameters["smrt_bin"]
    blasr = os.path.join(smrt_bin, "blasr")
    samtools = os.path.join( smrt_bin, "samtools")

    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , "aln_%s.sh" % (ctg_id))

    script = []
    script.append( "set -vex" )
    script.append( "trap 'touch {job_done}.exit' EXIT".format(job_done = job_done) )
    script.append( "cd %s" % wd )
    script.append( "hostname" )
    script.append( "date" )
    script.append( "cd {wd}".format(wd = wd) )
    script.append( "time {blasr} {read_fasta} {ref_fasta} -noSplitSubreads -clipping subread\
 -hitPolicy randombest -randomSeed 42 -bestn 1 -minPctIdentity 70.0\
 -minMatch 12  -nproc 24 -sam -out tmp_aln.sam".format(blasr = blasr,
                                                       read_fasta = read_fasta, 
                                                       ref_fasta = ref_fasta) )

    script.append( "{samtools} view -bS tmp_aln.sam | {samtools} sort - {ctg_id}_sorted".format( samtools = samtools, ctg_id = ctg_id) ) 
    script.append( "{samtools} index {ctg_id}_sorted.bam".format( samtools = samtools, ctg_id = ctg_id) )
    script.append( "rm tmp_aln.sam" )
    script.append( "date" )
    script.append( "touch {job_done}".format(job_done = job_done) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script) + '\n')

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_blasr_aln
    #run_script(job_data, job_type = config["job_type"])
    run_script(job_data, job_type = "SGE")
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])

def task_htig_asm(self):

    ref_fasta = fn(self.ref_fasta)
    aln_bam = fn(self.aln_bam)

    job_done = fn(self.job_done)

    job_uid = self.parameters["job_uid"]
    wd = self.parameters["wd"]
    sge_hasm = self.parameters["sge_hasm"]
    ctg_id = self.parameters["ctg_id"]

    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , "hasm_%s.sh" % (ctg_id))

    script = []

    script.append( "set -vex" )
    script.append( "trap 'touch {job_done}.exit' EXIT".format(job_done = job_done) )
    script.append( "cd %s" % wd )
    script.append( "hostname" )
    script.append( "date" )
    script.append( "cd {wd}".format(wd = wd) )
    script.append( "fc_phasing.py --bam {aln_bam} --fasta {ref_fasta} --ctg_id {ctg_id} --base_dir ../".format( aln_bam = aln_bam,
                                                                                                                ref_fasta = ref_fasta,
                                                                                                                ctg_id = ctg_id ))
    script.append( "fc_phasing_readmap.py --ctg_id {ctg_id} --read_map_dir ../../2-asm-falcon/read_maps --phased_reads phased_reads".format(ctg_id = ctg_id) )
    script.append( "fc_ovlp_filter_with_phase.py --fofn ../../2-asm-falcon/las.fofn\
            --max_diff 120 --max_cov 120 --min_cov 1 --n_core 12 --min_len 2500\
            --db ../../1-preads_ovl/preads.db  --rid_phase_map ./rid_to_phase > preads.p_ovl") #TODO: make it configurable
    script.append( "fc_phased_ovlp_to_graph.py preads.p_ovl --min_len 2500 > fc.log" )
    script.append( "fc_graphs_to_h_tigs.py --fc_asm_path ../../2-asm-falcon/ --fc_hasm_path ./ --ctg_id {ctg_id}\
            --rid_phase_map ./rid_to_phase --fasta ../../1-preads_ovl/preads4falcon.fasta".format(ctg_id = ctg_id))

    script.append( "fc_dedup_h_tigs.py" )
    #script.append( "cd ../../" )
    #script.append( "fc_track_reads_htigs.py {ctg_id}".format(ctg_id=ctg_id) )
    script.append( "date" )
    script.append( "touch {job_done}".format(job_done = job_done) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script) + '\n')

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_hasm
    #run_script(job_data, job_type = config["job_type"])
    run_script(job_data, job_type = "SGE")
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])



if __name__ == "__main__":
    global fc_run_logger
    fc_run_logger = support.setup_logger(None)
    ctg_ids = []
    with open("./3-unzip/reads/ctg_list") as f:
        for row in f:
            row = row.strip()
            ctg_ids.append( row )

    concurrent_jobs = 64
    PypeThreadWorkflow.setNumThreadAllowed(concurrent_jobs, concurrent_jobs)
    wf = PypeThreadWorkflow()

    sge_blasr_aln = " -pe smp 24 -q bigmem " 
    ctg_list_file = makePypeLocalFile("./3-unzip/reads/ctg_list")

    aln1_outs = {}
    for ctg_id in ctg_ids:
        # inputs
        ref_fasta = makePypeLocalFile("./3-unzip/reads/{ctg_id}_ref.fa".format(ctg_id = ctg_id))
        read_fasta = makePypeLocalFile("./3-unzip/reads/{ctg_id}_reads.fa".format(ctg_id = ctg_id))
        
        # outputs
        wd = os.path.join( os.getcwd(),  "./3-unzip/{ctg_id}/".format( ctg_id = ctg_id ) )
        mkdir(wd)
        ctg_aln_out = makePypeLocalFile( os.path.join( wd, "{ctg_id}_sorted.bam".format( ctg_id = ctg_id ) ) )
        job_done = makePypeLocalFile( os.path.join( wd, "aln_{ctg_id}_done".format( ctg_id = ctg_id ) ) )
    
        parameters = {"job_uid":"aln-"+ctg_id, "wd": wd, "sge_blasr_aln":sge_blasr_aln, "ctg_id": ctg_id,
                "smrt_bin":"/mnt/secondary/builds/full/3.0.0/prod/smrtanalysis_3.0.0.153854/smrtcmds/bin/"} 
        make_blasr_task = PypeTask(inputs = {"ref_fasta": ref_fasta, "read_fasta": read_fasta},
                                   outputs = {"ctg_aln_out": ctg_aln_out, "job_done": job_done},
                                   parameters = parameters,
                                   TaskType = PypeThreadTaskBase,
                                   URL = "task://localhost/aln1_{ctg_id}".format( ctg_id = ctg_id ) )
        blasr_task = make_blasr_task(task_run_blasr)
        aln1_outs[ctg_id] = (ctg_aln_out, job_done)
        wf.addTask(blasr_task)

        sge_hasm = " -pe smp 12 -q bigmem"
        job_done = makePypeLocalFile( os.path.join( wd, "hasm_{ctg_id}_done".format( ctg_id = ctg_id ) ) )
        parameters = {"job_uid":"ha-"+ctg_id, "wd": wd, "sge_hasm":sge_hasm, "ctg_id": ctg_id} 
        make_hasm_task = PypeTask(inputs = {"ref_fasta": ref_fasta, "aln_bam":ctg_aln_out},
                                   outputs = {"job_done": job_done},
                                   parameters = parameters,
                                   TaskType = PypeThreadTaskBase,
                                   URL = "task://localhost/h_{ctg_id}".format( ctg_id = ctg_id ) )
        hasm_task = make_hasm_task(task_htig_asm)
        wf.addTask(hasm_task)
    #print aln1_outs
    wf.refreshTargets()
        

