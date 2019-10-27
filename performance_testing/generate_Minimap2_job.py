# Given an input size, create a Minimap2 job script and launch
import argparse
import os
from string import Template

def get_options():
    """ Fetches the arguments for the program """

    program_desc = ("Given an input size, create a Minimap2 job script and "
                    " launch it.")
    parser = argparse.ArgumentParser(description=program_desc)

    parser.add_argument("--dir", "-d", dest = "in_dir", metavar='FILE',
        help = "Directory containing input sam files", type = str)
    parser.add_argument('--input', "-i", dest = 'input_size', type = str,
        help='Optional: comma-delimited input sizes to run', default = "all")
    parser.add_argument("--logs", "-l", dest = "logdir",
                        help = "Output directory for log files", type = str)
    parser.add_argument("--outdir", "-o", dest = "outdir", type = str,
                        help = "Output directory for all other run files")

    args = parser.parse_args()
    return args

def fill_template_script(n_reads, job_name, log_out, run_dir, fastq_file,
                        email = "dwyman@uci.edu", cores = 16,
                        minimap_path = "~/minimap2-2.15",
                        ref_genome = "/pub/dwyman/TALON-paper-2019/refs/hg38/hg38.fa"):

    template = Template(\
'''#!/bin/bash
#$$ -q sam128
#$$ -pe one-node-mpi $cores
#$$ -R y
#$$ -N $job_name
#$$ -M $email
#$$ -m ea
#$$ -o $log_out/
#$$ -j y

set -e
mkdir -p $run_dir
cd $run_dir

time $minimap_path/minimap2 \
    -t $cores \
    -ax splice -uf --secondary=no -C5 \
    $ref_genome \
    $fastq_file \
    > $run_dir/mapped.sam

source deactivate
qstat -j $$JOB_ID
echo input_size:$n_reads''')

    return template.substitute({'cores':cores,
                            'job_name':job_name,
                            'email':email,
                            'log_out':log_out,
                            'run_dir':run_dir,
                            'fastq_file':fastq_file,
                            'n_reads':n_reads,
                            'ref_genome':ref_genome,
                            'minimap_path': minimap_path})

def human_format(num):
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    num = int(num)
    # add more suffixes if you need them
    return '%d%s' % (num, ['', 'K', 'M', 'G', 'T', 'P'][magnitude])

def construct_job_name(n_reads, n_cores=16):
    """ Construct a name to be used in the job script based on the number of
        reads in the input as well as the number of cores """

    return "_".join(["Minimap", human_format(n_reads), str(n_cores)])

def construct_run_dir(n_reads, outdir):
    """ Create a directory for the run output based on the job name """

    run_dir = outdir + "/" + str(n_reads) + "/Minimap2"
    os.system("mkdir -p " + run_dir)
    return run_dir

def main():
    options = get_options()
    n_reads = int(options.input_size)
    fastq_file = options.in_dir + "/" + str(n_reads) + ".fq"
    job_name = construct_job_name(n_reads)
    run_dir = construct_run_dir(n_reads, options.outdir)
    log_out = "/".join([options.logdir, "Minimap2"]) 

    curr_script = run_dir + "/run_Minimap2.sh"
    with open(curr_script, 'w') as o:
        script_text = fill_template_script(n_reads, job_name, log_out, 
                                               run_dir, fastq_file)
        o.write(script_text)
    
    os.system("qsub %s" % (curr_script))

if __name__ == '__main__':
    main()
