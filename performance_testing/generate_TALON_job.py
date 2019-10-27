# Given an input size, create a TALON job script and launch
import argparse
import os
from string import Template

def get_options():
    """ Fetches the arguments for the program """

    program_desc = ("Given an input size, create a TALON job script and "
                    "launch it.")
    parser = argparse.ArgumentParser(description=program_desc)

    parser.add_argument('--input', "-i", dest = 'input_size', type = str,
        help='Optional: comma-delimited input sizes to run', default = "all")
    parser.add_argument("--logs", "-l", dest = "logdir",
                        help = "Output directory for log files", type = str)
    parser.add_argument("--outdir", "-o", dest = "outdir", type = str,
                        help = "Output directory for all other run files")

    args = parser.parse_args()
    return args

def make_config(read_count, run_dir, sam_file):
    """ Create a config file for run """
    config_file = run_dir + "/config.csv"
    dataset = "test"
    metadata = "test"
    platform = "PacBio"

    with open(config_file, 'w') as o:
        o.write(",".join([dataset, metadata, platform, sam_file]) + "\n")

    return config_file

def fill_template_script(n_reads, job_name, log_out, run_dir, config_file,
                         email = "dwyman@uci.edu", cores = 16,
                         db_path = "/pub/dwyman/TALON-paper-2019/refs/TALON/unmodified_full_gencode_v29_2019-06-19.db", conda_env = "python3.6"):

    # Make a config file

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
source activate $conda_env
module load samtools
mkdir -p $run_dir
mkdir -p /scratch/$job_name
cd $run_dir
cp $db_path talon.db

time talon --f $config_file \\
           --db talon.db \\
           --build hg38 \\
           --cov 0.9 \\
           --identity 0 \\
           --o $run_dir/talon

source deactivate
qstat -j $$JOB_ID
echo input_size:$n_reads''')

    return template.substitute({'cores':cores,
                            'job_name':job_name,
                            'config_file':config_file,
                            'email':email,
                            'log_out':log_out,
                            'run_dir':run_dir,
                            'n_reads':n_reads,
                            'db_path': db_path,
                            'conda_env':conda_env})

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

    return "_".join(["TALON", human_format(n_reads), str(n_cores)])

def construct_run_dir(n_reads, outdir):
    """ Create a directory for the run output based on the job name """

    run_dir = outdir + "/" + str(n_reads) + "/TALON"
    os.system("mkdir -p " + run_dir)
    return run_dir

def main():
    options = get_options()
    n_reads = int(options.input_size)
    sam_file = options.outdir + "/" + str(n_reads) + "/TC/TC_clean.sam"
    job_name = construct_job_name(n_reads)
    run_dir = construct_run_dir(n_reads, options.outdir)
    log_out = "/".join([options.logdir, "TALON"]) 

    curr_script = run_dir + "/run_TALON.sh"
    with open(curr_script, 'w') as o:
        config_file = make_config(n_reads, run_dir, sam_file)
        script_text = fill_template_script(n_reads, job_name, log_out, run_dir,                                            config_file)
        o.write(script_text)
    
    #os.system("qsub %s" % (curr_script))

if __name__ == '__main__':
    main()
