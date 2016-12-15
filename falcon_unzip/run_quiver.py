from falcon_kit import run_support as support
from pypeflow.simple_pwatcher_bridge import (
        PypeLocalFile, makePypeLocalFile, fn,
        PypeTask,
        PypeProcWatcherWorkflow, MyFakePypeThreadTaskBase)
from falcon_kit.FastaReader import FastaReader
import glob
import json
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


def task_run_quiver(self):
    ref_fasta = fn(self.ref_fasta)
    read_sam = fn(self.read_sam)

    cns_fasta = fn(self.cns_fasta)
    cns_fastq = fn(self.cns_fastq)
    job_done = fn(self.job_done)

    job_uid = self.parameters['job_uid']
    ctg_id = self.parameters['ctg_id']

    smrt_bin = self.parameters['smrt_bin']
    samtools = os.path.join(smrt_bin, 'samtools')
    pbalign = os.path.join(smrt_bin, 'pbalign')
    makePbi = os.path.join(smrt_bin, 'makePbi')
    variantCaller = os.path.join(smrt_bin, 'variantCaller')

    script_fn = 'cns_%s.sh' % (ctg_id)

    script = """\
set -vex
trap 'touch {job_done}.exit' EXIT
hostname
date

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

def rm(f):
    system('rm -f {}'.format(f))
def touch(f):
    system('touch {}'.format(f))
def task_cns_zcat(self):
    gathered_p_ctg_fn = fn(self.gathered_p_ctg)
    gathered_h_ctg_fn = fn(self.gathered_h_ctg)

    cns_p_ctg_fasta_fn = fn(self.cns_p_ctg_fasta)
    cns_p_ctg_fastq_fn = fn(self.cns_p_ctg_fastq)
    cns_h_ctg_fasta_fn = fn(self.cns_h_ctg_fasta)
    cns_h_ctg_fastq_fn = fn(self.cns_h_ctg_fastq)
    job_done_fn = fn(self.job_done)

    rm(cns_p_ctg_fasta_fn)
    touch(cns_p_ctg_fasta_fn)
    rm(cns_p_ctg_fastq_fn)
    touch(cns_p_ctg_fastq_fn)
    with open(gathered_p_ctg_fn) as ifs:
        for line in ifs:
            cns_fasta_fn, cns_fastq_fn = line.split()
            system('zcat {cns_fasta_fn} >> {cns_p_ctg_fasta_fn}'.format(**locals()))
            system('zcat {cns_fastq_fn} >> {cns_p_ctg_fastq_fn}'.format(**locals()))
    with open(gathered_p_ctg_fn) as ifs:
        for line in ifs:
            cns_fasta_fn, cns_fastq_fn = line.split()
            rm(cns_fasta_fn)
            rm(cns_fasta_fn)
    rm(cns_h_ctg_fasta_fn)
    touch(cns_h_ctg_fasta_fn)
    rm(cns_h_ctg_fastq_fn)
    touch(cns_h_ctg_fastq_fn)
    with open(gathered_h_ctg_fn) as ifs:
        for line in ifs:
            cns_fasta_fn, cns_fastq_fn = line.split()
            system('zcat {cns_fasta_fn} >> {cns_h_ctg_fasta_fn}'.format(**locals()))
            system('zcat {cns_fastq_fn} >> {cns_h_ctg_fastq_fn}'.format(**locals()))
    with open(gathered_h_ctg_fn) as ifs:
        for line in ifs:
            cns_fasta_fn, cns_fastq_fn = line.split()
            rm(cns_fasta_fn)
            rm(cns_fasta_fn)

    touch(job_done_fn)

def task_scatter_quiver(self):
    p_ctg_fn = fn(self.p_ctg_fa)
    h_ctg_fn = fn(self.h_ctg_fa)
    out_json = fn(self.scattered_quiver_json)
    config = self.parameters['config']

    ref_seq_data = {}

    # I think this will crash if the file is empty. Maybe that is ok.
    p_ctg_fa = FastaReader(p_ctg_fn)
    ctg_types = {}
    for r in p_ctg_fa:
        rid = r.name.split()[0]
        ref_seq_data[rid] = r.sequence
        ctg_types[rid] = 'p'


    # I think this will crash if the file is empty. Maybe that is ok.
    h_ctg_fa = FastaReader(h_ctg_fn)
    for r in h_ctg_fa:
        rid = r.name.split()[0]
        ref_seq_data[rid] = r.sequence
        ctg_types[rid] = 'h'

    ctg_ids = sorted(ref_seq_data.keys())
    #p_ctg_out=[]
    #h_ctg_out=[]
    #job_done_plfs = {}
    jobs = []
    for ctg_id in ctg_ids:
        sequence = ref_seq_data[ctg_id]
        m_ctg_id = ctg_id.split('-')[0]
        wd = os.path.join(os.getcwd(), m_ctg_id)
        ref_fasta = os.path.join(wd, '{ctg_id}_ref.fa'.format(ctg_id = ctg_id))
        read_sam = os.path.join(os.getcwd(), '../reads/{ctg_id}.sam'.format(ctg_id = ctg_id))
        #cns_fasta = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fasta.gz'.format(ctg_id = ctg_id)))
        #cns_fastq = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fastq.gz'.format(ctg_id = ctg_id)))
        #job_done = makePypeLocalFile(os.path.join(wd, '{ctg_id}_quiver_done'.format(ctg_id = ctg_id)))

        if os.path.exists(read_sam): # TODO(CD): Ask Jason what we should do if missing SAM. And what about network latency?
            #if ctg_types[ctg_id] == 'p':
            #    p_ctg_out.append( (cns_fasta, cns_fastq) )
            #if ctg_types[ctg_id] == 'h':
            #    h_ctg_out.append( (cns_fasta, cns_fastq) )
            mkdir(wd)
            if not os.path.exists(ref_fasta):
                with open(ref_fasta,'w') as f:
                    print >>f, '>'+ctg_id
                    print >>f, sequence
            #parameters = {'job_uid':'q-'+ctg_id, 'wd': wd, 'config':config, 'ctg_id': ctg_id }
            #make_quiver_task = PypeTask(inputs = {'ref_fasta': ref_fasta, 'read_sam': read_sam},
            #                           outputs = {'cns_fasta': cns_fasta, 'cns_fastq': cns_fastq, 'job_done': job_done},
            #                           parameters = parameters,
            #                           )
            #quiver_task = make_quiver_task(task_run_quiver)
            #wf.addTask(quiver_task)
            #job_done_plfs['{}'.format(ctg_id)] = job_done
            new_job = {}
            new_job['ctg_id'] = ctg_id
            new_job['ctg_types'] = ctg_types
            new_job['smrt_bin'] = config['smrt_bin']
            new_job['sge_option'] = config['sge_quiver']
            new_job['ref_fasta'] = ref_fasta
            new_job['read_sam'] = read_sam
            jobs.append(new_job)
    open(out_json, 'w').write(json.dumps(jobs))

def create_quiver_jobs(wf, scattered_quiver_plf):
    scattered_quiver_fn = fn(scattered_quiver_plf)
    jobs = json.loads(open(scattered_quiver_fn).read())
    #ctg_ids = sorted(jobs['ref_seq_data'])
    p_ctg_out=[]
    h_ctg_out=[]
    job_done_plfs = {}
    for job in jobs:
        ctg_id = job['ctg_id']
        ctg_types = job['ctg_types']
        smrt_bin = job['smrt_bin']
        sge_option = job['sge_option']
        ref_fasta = makePypeLocalFile(job['ref_fasta'])
        read_sam = makePypeLocalFile(job['read_sam'])
        m_ctg_id = ctg_id.split('-')[0]
        wd = os.path.join(os.getcwd(), './4-quiver/', m_ctg_id)
        #ref_fasta = makePypeLocalFile(os.path.join(wd, '{ctg_id}_ref.fa'.format(ctg_id = ctg_id)))
        #read_sam = makePypeLocalFile(os.path.join(os.getcwd(), './4-quiver/reads/' '{ctg_id}.sam'.format(ctg_id = ctg_id)))
        cns_fasta = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fasta.gz'.format(ctg_id = ctg_id)))
        cns_fastq = makePypeLocalFile(os.path.join(wd, 'cns-{ctg_id}.fastq.gz'.format(ctg_id = ctg_id)))
        job_done = makePypeLocalFile(os.path.join(wd, '{ctg_id}_quiver_done'.format(ctg_id = ctg_id)))

        if os.path.exists(fn(read_sam)): # TODO(CD): Ask Jason what we should do if missing SAM.
            if ctg_types[ctg_id] == 'p':
                p_ctg_out.append( (fn(cns_fasta), fn(cns_fastq)) )
            elif ctg_types[ctg_id] == 'h':
                h_ctg_out.append( (fn(cns_fasta), fn(cns_fastq)) )
            else:
                LOG.warning('Type is {!r}, not "p" or "h". Why are we running Quiver?'.format(ctg_types[ctg_id]))
            parameters = {
                    'job_uid':'q-'+ctg_id,
                    'ctg_id': ctg_id,
                    'smrt_bin': smrt_bin,
                    'sge_option': sge_option,
            }
            make_quiver_task = PypeTask(inputs = {'ref_fasta': ref_fasta, 'read_sam': read_sam,
                                         'scattered_quiver': scattered_quiver_plf,
                                       },
                                       outputs = {'cns_fasta': cns_fasta, 'cns_fastq': cns_fastq, 'job_done': job_done},
                                       parameters = parameters,
            )
            quiver_task = make_quiver_task(task_run_quiver)
            wf.addTask(quiver_task)
            job_done_plfs['{}'.format(ctg_id)] = job_done
    return p_ctg_out, h_ctg_out, job_done_plfs

def task_gather_quiver(self):
    """We wrote the "gathered" files during task construction.
    """
    job_done_fn = fn(self.job_done)
    touch(job_done_fn)


def main(argv=sys.argv):
    global LOG
    LOG = support.setup_logger(None)


    if len(sys.argv) < 2:
        print>>sys.stderr, 'you need to provide a configuration file to specific a couple cluster running environment'
        sys.exit(1)

    config_fn = sys.argv[1]
    config_absbasedir = os.path.dirname(os.path.abspath(config_fn))

    config = ConfigParser.ConfigParser()
    config.read(config_fn)


    job_type = 'SGE'
    if config.has_option('General', 'job_type'):
        job_type = config.get('General', 'job_type')

    sge_track_reads = ' -pe smp 12 -q bigmem'
    if config.has_option('Unzip', 'sge_track_reads'):
        sge_track_reads = config.get('Unzip', 'sge_track_reads')

    sge_quiver = ' -pe smp 24 -q bigmem '
    if config.has_option('Unzip', 'sge_quiver'):
        sge_quiver = config.get('Unzip', 'sge_quiver')

    smrt_bin = '/mnt/secondary/builds/full/3.0.0/prod/smrtanalysis_3.0.0.153854/smrtcmds/bin/'
    if config.has_option('Unzip', 'smrt_bin'):
        smrt_bin = config.get('Unzip', 'smrt_bin')

    input_bam_fofn = 'input_bam.fofn'
    if config.has_option('Unzip', 'input_bam_fofn'):
        input_bam_fofn = config.get('Unzip', 'input_bam_fofn')
    if not os.path.isabs(input_bam_fofn):
        input_bam_fofn = os.path.join(config_absbasedir, input_bam_fofn)


    quiver_concurrent_jobs = 8
    if config.has_option('Unzip', 'quiver_concurrent_jobs'):
        quiver_concurrent_jobs = config.getint('Unzip', 'quiver_concurrent_jobs')

    config = {'job_type': job_type,
              'sge_quiver': sge_quiver,
              'sge_track_reads': sge_track_reads,
              'input_bam_fofn': input_bam_fofn,
              'smrt_bin': smrt_bin}
    LOG.info('config={}'.format(pprint.pformat(config)))

    #support.job_type = 'SGE' #tmp hack until we have a configuration parser


    wf = PypeProcWatcherWorkflow(
            max_jobs=quiver_concurrent_jobs,
            job_type=config['job_type'],
            job_queue=config.get('job_queue'),
            sge_option=config.get('sge_option'),
            watcher_type=config.get('pwatcher_type'),
            #watcher_directory=config.get('pwatcher_directory', 'mypwatcher'),
            use_tmpdir=config.get('use_tmpdir'),
    )

    abscwd = os.path.abspath('.')
    parameters = {'wd': os.path.join(abscwd, '4-quiver', 'track_reads_h'), 'config': config,
            'sge_option': config['sge_track_reads'],
    }
    hasm_done_plf = makePypeLocalFile('./3-unzip/1-hasm/hasm_done') # by convention
    track_reads_h_done_plf = makePypeLocalFile(os.path.join(parameters['wd'], 'track_reads_h_done'))
    make_track_reads_task = PypeTask(inputs = {'hasm_done': hasm_done_plf},
                                     outputs = {'job_done': track_reads_h_done_plf},
                                     parameters = parameters,
    )
    track_reads_task = make_track_reads_task(task_track_reads)

    wf.addTask(track_reads_task)

    scattered_quiver_plf = makePypeLocalFile('4-quiver/quiver_scatter/scattered.json')
    parameters = {
            'config': config,
    }
    make_task = PypeTask(
            inputs = {
                'p_ctg_fa': makePypeLocalFile('3-unzip/all_p_ctg.fa'),
                'h_ctg_fa': makePypeLocalFile('3-unzip/all_h_ctg.fa'),
                'track_reads_h_done': track_reads_h_done_plf,
            },
            outputs = {
                'scattered_quiver_json': scattered_quiver_plf,
            },
            parameters = parameters,
    )
    wf.addTask(make_task(task_scatter_quiver))
    wf.refreshTargets()

    p_ctg_out, h_ctg_out, job_done_plfs = create_quiver_jobs(wf, scattered_quiver_plf)

    gathered_p_ctg_plf = makePypeLocalFile('4-quiver/cns_gather/p_ctg.txt')
    gathered_h_ctg_plf = makePypeLocalFile('4-quiver/cns_gather/h_ctg.txt')
    gather_done_plf = makePypeLocalFile('4-quiver/cns_gather/job_done')
    mkdir('4-quiver/cns_gather')
    with open(fn(gathered_p_ctg_plf), 'w') as ifs:
        for cns_fasta_fn, cns_fastq_fn in sorted(p_ctg_out):
            ifs.write('{} {}\n'.format(cns_fasta_fn, cns_fastq_fn))
    with open(fn(gathered_h_ctg_plf), 'w') as ifs:
        for cns_fasta_fn, cns_fastq_fn in sorted(h_ctg_out):
            ifs.write('{} {}\n'.format(cns_fasta_fn, cns_fastq_fn))

    make_task = PypeTask(
            inputs = job_done_plfs,
            outputs = {
                'job_done': gather_done_plf,
            },
            parameters = {},
    )
    wf.addTask(make_task(task_gather_quiver))
    wf.refreshTargets()

    cns_p_ctg_fasta_plf = makePypeLocalFile('4-quiver/cns_output/cns_p_ctg.fasta')
    cns_p_ctg_fastq_plf = makePypeLocalFile('4-quiver/cns_output/cns_p_ctg.fastq')
    cns_h_ctg_fasta_plf = makePypeLocalFile('4-quiver/cns_output/cns_h_ctg.fasta')
    cns_h_ctg_fastq_plf = makePypeLocalFile('4-quiver/cns_output/cns_h_ctg.fastq')
    zcat_done_plf = makePypeLocalFile('4-quiver/cns_output/job_done')
    make_task = PypeTask(
            inputs = {
                'gathered_p_ctg': gathered_p_ctg_plf,
                'gathered_h_ctg': gathered_h_ctg_plf,
                'gather_done': gather_done_plf,
            },
            outputs = {
                'cns_p_ctg_fasta': cns_p_ctg_fasta_plf,
                'cns_p_ctg_fastq': cns_p_ctg_fastq_plf,
                'cns_h_ctg_fasta': cns_h_ctg_fasta_plf,
                'cns_h_ctg_fastq': cns_h_ctg_fastq_plf,
                'job_done': zcat_done_plf,
            },
    )
    wf.addTask(make_task(task_cns_zcat))

    wf.refreshTargets()
