import argparse
from glob import glob
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import ScalarFormatter

def get_options():
    """ Fetches the arguments for the program """

    program_desc = ("Collects runtime, memory usage, and input size info from "
                    "TranscriptClean log files and plots the results.")

    parser = argparse.ArgumentParser(description=program_desc)

    parser.add_argument("--logs", "-l", dest = "logdir",
                        help = "Output directory for log files", type = str)
    parser.add_argument("--outprefix", "-o", dest = "outprefix", type = str,
                        help = "Output prefix for plots", default = "performance")

    args = parser.parse_args()
    return args

def process_logs_of_type(log_dir, program_name):
    """ Iterate through logs for a particular program and store performance
        metrics in a dataframe"""

    input_sizes = []
    runtimes = []
    memory_use = []
    n_cores = []

    log_dir = log_dir + "/" + program_name + "/"
    for log in glob("/".join([log_dir, "*.o*"])):
        try:
            input_size, runtime, memory, cores = parse_log(log)
            input_sizes.append(input_size)
            runtimes.append(runtime)
            memory_use.append(memory)
            n_cores.append(cores)
        except:
            print("Problem parsing file %s. Skipping..." % (log))

    # Now format the data for plotting
    df = pd.DataFrame({'n_reads': input_sizes,
                       'runtime': runtimes,
                       'memory': memory_use,
                       'n_cores': n_cores})
    df['program'] = program_name
    return df 

def parse_log(log):
    """Parse the provided log file to extract the input size ( reads), runtime,
    memory usage, and n_cores of the job."""

    with open(log, 'r') as f:
        log_contents = f.readlines()

        # Get number of reads in input (units of a thousand)
        input_size = int(log_contents[-1].split("input_size:")[-1].strip())/1000000

        # Get runtime in minutes
        runtime_line = [i for i in log_contents if i.startswith('real')][0]
        runtime_str = runtime_line.split()[-1]
        minutes = int(runtime_str.split("m")[0])
        seconds = float((runtime_str.split("m")[-1]).split("s")[0])
        runtime_m = round(minutes + seconds/60., 0)

        # Get memory usage
        maxvmem_line = [i for i in log_contents if i.startswith('usage')][0]
        maxvmem_str = maxvmem_line.split("maxvmem=")[-1]
        if "G" in maxvmem_str:
            maxvmem = round(float(maxvmem_str.split("G")[0]))
        elif "M" in maxvmem_str:
            maxvmem = float(maxvmem_str.split("M")[0])*0.001
        else:
            print("Memory not in units of G or M")

        # Get number of cores
        cores_line = [i for i in log_contents if i.startswith('parallel environment')][0]
        cores = int(cores_line.split("one-node-mpi range: ")[-1].strip())

    return input_size, runtime_m, maxvmem, cores

def plot_reads_vs_var(df, y_var, ylabel, fname):
    """ Plot read counts against another variable """

    xlabel = "Mapped read count at pipeline start (millions)"

    style = dict(size=16, color='black')
    sb.set_context("paper", font_scale=1.25)
    ax = sb.pointplot(data = df, x='n_reads', y=y_var, scale = 0.75)
    ax.set(xlabel=xlabel, ylabel=ylabel)
    #ax.set_yscale('log', basey=2)
    ax.grid(b=True, which='major')
    ax.set_axisbelow(True)
    ax.minorticks_on()

    plot = ax.get_figure()
    plot.savefig(fname, dpi = 600)
    plt.close(plot)

def main():

    options = get_options()
    log_dir = options.logdir
    outprefix = options.outprefix

    # Iterate over log files and process each to extract input and runtime
    minimap = process_logs_of_type(log_dir, "Minimap2")
    TC = process_logs_of_type(log_dir, "TC")
    TALON = process_logs_of_type(log_dir, "TALON")

    all_data = pd.concat([minimap, TC, TALON])
    
    total_runtimes = all_data.groupby("n_reads").agg({'runtime': 'sum'})
    plot_reads_vs_var(minimap, "runtime", "Runtime (min)", outprefix + "_minimap2_runtime.png")

if __name__ == '__main__':
    main()
