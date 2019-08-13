from glob import glob
import seaborn as sb
import numpy as np
import matplotlib.pyplot as plt

def get_read_count(curr_dir):
    """ Return in units of millions """
    rc_file = curr_dir + "read_count.txt"
    try:
        with open(rc_file) as f:
            first_line = f.readline()
            rc = round(int(first_line.strip())/1000000, 3)

    except:
        raise NameError("Problem processing file with name '%s'" % (rc_file))
        exit()

    return rc

def get_runtime_and_maxRAM(dataset):
    """ Return runtime in hours and max RAM in GB """

    curr_dir = dataset + "/"
    log = {}
    logfile = curr_dir + dataset + ".log"
    with open(logfile) as f:
        for line in f:
            if line.startswith("="):
                continue
            (key, val) = line.split()[0:2]
            log[key] = val
    
    runtime_s = float((log["ru_utime"]).split("s")[0])
    max_RAM = float(log["maxvmem"].split("GB")[0])

    # Convert to hours
    runtime = round(runtime_s/(60*60),3)
    return runtime, max_RAM

def plot_reads_vs_var(read_counts, y, ylabel, fname):
    """ Plot read counts against another variable """

    xlabel = "Read count (millions)"

    style = dict(size=16, color='black')
    sb.set_context("paper", font_scale=1.5)
    ax = sb.regplot(x = read_counts, y = y)
    ax.set(xlabel=xlabel, ylabel=ylabel)

    plot = ax.get_figure()
    #corrfunc(x, y)
    plot.savefig(fname, dpi = 600)
    plt.close(plot)    


def main():
    # Glob all subdirectories- each one contains a run
    datasets = [x.split("/")[1] for x in glob("./*/")]

    # For each subdir, get the read count and parse the runtime and max RAM
    # out of the log file
    read_counts = []
    runtimes = []
    max_RAM = []

    for dataset in datasets:
        curr_dir = dataset + "/"
    
        # Get read count
        rc = get_read_count(curr_dir)
        read_counts.append(rc)

        # Get runtime and max RAM
        rt, mr = get_runtime_and_maxRAM(dataset)
        runtimes.append(rt)
        max_RAM.append(mr)

    print(datasets)
    print(read_counts)
    print(runtimes)
    print(max_RAM)

    # Plot the read count vs. the runtime
    plot_reads_vs_var(read_counts, runtimes, "Runtime (hours)", "rc_vs_runtime.png")

    # Plot the read count vs. the maximum RAM used during the run (maxvmem)
    plot_reads_vs_var(read_counts, max_RAM, "Maximum RAM usage (GB)", "rc_vs_maxRAM.png")

if __name__ == '__main__':
    main()


