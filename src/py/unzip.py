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
#fc_run_logger = None

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

def task_track_reads(self):

    job_done = fn(self.job_done)
    wd = self.parameters["wd"]
    config = self.parameters["config"]
    sge_track_reads = config["sge_track_reads"]
    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , "track_reads.sh")
    
    script = []
    script.append( "set -vex" )
    script.append( "trap 'touch {job_done}.exit' EXIT".format(job_done = job_done) )
    script.append( "cd %s" % wd )
    script.append( "hostname" )
    script.append( "date" )
    script.append( "cd {wd}".format(wd = wd) )
    script.append( "python -m falcon_kit.mains.get_read_ctg_map" )
    script.append( "python -m falcon_kit.mains.rr_ctg_track" )
    script.append( "python -m falcon_kit.mains.pr_ctg_track" )
    script.append( "mkdir -p 3-unzip/reads/" )
    script.append( "python -m falcon_kit.mains.fetch_reads" )
    script.append( "date" )
    script.append( "touch {job_done}".format(job_done = job_done) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script) + '\n')

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_track_reads
    run_script(job_data, job_type = config["job_type"])
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])


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
    run_script(job_data, job_type = config["job_type"])
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])

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
    script.append( "fc_phasing_readmap.py --ctg_id {ctg_id} --read_map_dir ../../../2-asm-falcon/read_maps --phased_reads phased_reads".format(ctg_id = ctg_id) )
    script.append( "date" )
    script.append( "touch {job_done}".format(job_done = job_done) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script) + '\n')

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_phasing
    run_script(job_data, job_type = job_type)
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])

def task_hasm(self):

    job_done = fn(self.job_done)
    config = self.parameters["config"]
    sge_hasm = config["sge_hasm"]

    wd = self.parameters["wd"]

    job_type = config["job_type"]

    script_dir = os.path.join( wd )
    script_fn =  os.path.join( script_dir , "hasm.sh")

    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date
cd {wd}

fc_ovlp_filter_with_phase.py --fofn ../../2-asm-falcon/las.fofn --max_diff 120 --max_cov 120 --min_cov 1 --n_core 12 --min_len 2500 --db ../../1-preads_ovl/preads.db --rid_phase_map ./rid_to_phase.all > preads.p_ovl
fc_phased_ovlp_to_graph.py preads.p_ovl --min_len 2500 > fc.log
ln -sf db2falcon/preads4falcon.fasta ../../1-preads_ovl/
fc_graphs_to_h_tigs.py --fc_asm_path ../../2-asm-falcon/ --fc_hasm_path ./ --ctg_id all --rid_phase_map ./rid_to_phase.all --fasta ../../1-preads_ovl/preads4falcon.fasta
""".format(**locals())
    more_script = \
"""
WD=$PWD
for f in `cat ../reads/ctg_list `;do cd $WD/$f; fc_dedup_h_tigs.py $f; done

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
"""
    script.append( more_script )  # a little bit hacky here, we should improve
    script.append( "date" )
    script.append( "touch {job_done}".format(job_done = job_done) )

    with open(script_fn,"w") as script_file:
        script_file.write("\n".join(script) + '\n')

    job_data = support.make_job_data(self.URL, script_fn)
    job_data["sge_option"] = sge_hasm
    run_script(job_data, job_type = job_type)
    wait_for_file(job_done, task=self, job_name=job_data['job_name'])


def unzip_all(config):
    unzip_concurrent_jobs = config["unzip_concurrent_jobs"]
    PypeThreadWorkflow.setNumThreadAllowed(unzip_concurrent_jobs, unzip_concurrent_jobs)
    wf = PypeThreadWorkflow()

    ctg_list_file = makePypeLocalFile("./3-unzip/reads/ctg_list")
    falcon_asm_done = makePypeLocalFile("./2-asm-falcon/falcon_asm_done")
    parameters = {"wd": os.path.abspath("."), "config": config}

    job_done = makePypeLocalFile( os.path.join( parameters["wd"], "track_reads_done" ) )
    make_track_reads_task = PypeTask(inputs = {"falcon_asm_done": falcon_asm_done},
                                     outputs = {"job_done": job_done, "ctg_list_file": ctg_list_file},
                                     parameters = parameters,
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
        mkdir(wd)
        ctg_aln_out = makePypeLocalFile( os.path.join( wd, "{ctg_id}_sorted.bam".format( ctg_id = ctg_id ) ) )
        job_done = makePypeLocalFile( os.path.join( wd, "aln_{ctg_id}_done".format( ctg_id = ctg_id ) ) )
    
        parameters = {"job_uid":"aln-"+ctg_id, "wd": wd, "config":config, "ctg_id": ctg_id} 
        make_blasr_task = PypeTask(inputs = {"ref_fasta": ref_fasta, "read_fasta": read_fasta},
                                   outputs = {"ctg_aln_out": ctg_aln_out, "job_done": job_done},
                                   parameters = parameters,
                                   TaskType = PypeThreadTaskBase,
                                   URL = "task://localhost/aln_{ctg_id}".format( ctg_id = ctg_id ) )
        blasr_task = make_blasr_task(task_run_blasr)
        aln1_outs[ctg_id] = (ctg_aln_out, job_done)
        wf.addTask(blasr_task)

        job_done = makePypeLocalFile( os.path.join( wd, "p_{ctg_id}_done".format( ctg_id = ctg_id ) ) )
        rid_to_phase_out = makePypeLocalFile( os.path.join( wd, "rid_to_phase.{ctg_id}".format( ctg_id = ctg_id ) ) )
        all_ctg_out[ "r2p.{ctg_id}".format( ctg_id = ctg_id ) ] = rid_to_phase_out

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
    mkdir(hasm_wd)
    rid_to_phase_all =  makePypeLocalFile( os.path.join(hasm_wd, "rid_to_phase.all") )
    @PypeTask(inputs = all_ctg_out, outputs = {"rid_to_phase_all": rid_to_phase_all},
                TaskType = PypeThreadTaskBase, URL = "task://localhost/rid_to_phase_all" )
    def get_rid_to_phase_all(self):
        rid_to_phase_all_fn = fn(self.rid_to_phase_all)
        inputs_fn = [ fn(f) for f in self.inputs.values() ]
        inputs_fn.sort()
        output = []
        for fname in inputs_fn:
            output.extend( open(fname).read() )

        out = open( rid_to_phase_all_fn, "w")
        out.write("".join(output))
        out.close()

    wf.addTask(get_rid_to_phase_all)

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
        

def main(argv=sys.argv):

    if len(argv) < 2:
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

    support.job_type = "SGE" #tmp hack until we have a configuration parser

    unzip_all(config)
