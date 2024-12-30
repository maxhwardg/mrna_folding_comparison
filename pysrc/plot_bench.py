import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import numpy as np


def main():
    # Parse the arguments
    parser = ap.ArgumentParser(description="Plot the benchmark data")
    # Take any number of paths to data files
    parser.add_argument(
        "--data",
        "-d",
        nargs="+",
        help="Path to the data file. Multiple files can be given."
        "If multiple files are given, they are all plotted"
        " together with high/low/median shown.",
        required=True,
    )
    parser.add_argument(
        "--scale", "-s", help="Scale of the plot (linear/log)", default="linear"
    )
    args = parser.parse_args()

    assert args.scale in ["linear", "log"], "Invalid scale. Choose from linear or log"
    plot_scale = args.scale

    # Check for the data files
    for path in args.data:
        try:
            file = open(path, "r")
        except FileNotFoundError:
            print(f"File not found: {path}")
            return

    file = open("../data/bench_runs/randomseq_benchmark_data", "r")
    raw_data = file.read()
    plot_scale = "linear"

    # Function to parse the raw data
    def parse_data(raw_data):
        data = {
            "aa_len": [],
            "lineardesign_time": [],
            "cdsfold_time": [],
            "derna_time": [],
            "mrnafold_time": [],
            "lineardesign_mem": [],
            "cdsfold_mem": [],
            "derna_mem": [],
            "mrnafold_mem": [],
        }
        for line in raw_data.split("\n"):
            if line.startswith("aa_len:"):
                data["aa_len"].append(int(line.split(":")[1].strip()))
            elif line.startswith("lineardesign time(s):"):
                data["lineardesign_time"].append(float(line.split(":")[1].strip()))
            elif line.startswith("lineardesign memory(bytes):"):
                data["lineardesign_mem"].append(int(line.split(":")[1].strip()))
            elif line.startswith("cdsfold time(s):"):
                data["cdsfold_time"].append(float(line.split(":")[1].strip()))
            elif line.startswith("cdsfold memory(bytes):"):
                data["cdsfold_mem"].append(int(line.split(":")[1].strip()))
            elif line.startswith("derna time(s):"):
                data["derna_time"].append(float(line.split(":")[1].strip()))
            elif line.startswith("derna memory(bytes):"):
                data["derna_mem"].append(int(line.split(":")[1].strip()))
            elif line.startswith("mrnafold time(s):"):
                data["mrnafold_time"].append(float(line.split(":")[1].strip()))
            elif line.startswith("mrnafold memory(bytes):"):
                data["mrnafold_mem"].append(int(line.split(":")[1].strip()))
            else:
                print(f"Skipping line: {line}")
                continue

        def convert_to_gb(lst):
            return [None if x is None else x / 1024**3 for x in lst]

        data["lineardesign_mem"] = convert_to_gb(data["lineardesign_mem"])
        data["cdsfold_mem"] = convert_to_gb(data["cdsfold_mem"])
        data["derna_mem"] = convert_to_gb(data["derna_mem"])
        data["mrnafold_mem"] = convert_to_gb(data["mrnafold_mem"])
        
        
        def convert_to_minutes(lst):
            return [None if x is None else x / 60 for x in lst]
        
        data["lineardesign_time"] = convert_to_minutes(data["lineardesign_time"])
        data["cdsfold_time"] = convert_to_minutes(data["cdsfold_time"])
        data["derna_time"] = convert_to_minutes(data["derna_time"])
        data["mrnafold_time"] = convert_to_minutes(data["mrnafold_time"])

        return data

    # Parse the data
    parsed_data = []
    for path in args.data:
        file = open(path, "r")
        raw_data = file.read()
        parsed_data.append(parse_data(raw_data))
        # Check the same set of aa_len values are used
        if parsed_data[0]["aa_len"] != parsed_data[-1]["aa_len"]:
            print("Mismatch in aa_len values between files. Exiting.")
            return
        print(f"Data parsed from {path}")
        
    df_parsed = pd.DataFrame(parsed_data)
    print(df_parsed)
    
    
    
    plot_labels = ["LinearDesign", "CDSfold", "DERNA", "mRNAfold"]
    markers = ["o", "s", "^", "d"]
    time_data_labels = ["lineardesign_time", "cdsfold_time", "derna_time", "mrnafold_time"]
    mem_data_labels = ["lineardesign_mem", "cdsfold_mem", "derna_mem", "mrnafold_mem"]
    
    mins, maxs = df_parsed.min(axis=0), df_parsed.max(axis=0)
    def get_min_max_med(label):
        mat = []
        for row in df_parsed[label]:
            rvals = []
            for val in row:
                if val is not None:
                    rvals.append(val)
            mat.append(rvals)
        # transpose the matrix
        mat = list(map(list, zip(*mat)))
        meds = []
        mins = []
        maxs = []
        for row in mat:
            row.sort()
            meds.append(row[len(row) // 2])
            mins.append(row[0])
            maxs.append(row[-1])
        return mins, maxs, meds
            
    
    def make_plot(data_labels):
        plt.figure(figsize=(12, 8))
        for i, plot_label in enumerate(plot_labels):
            marker = markers[i]
            data_label = data_labels[i]
            mins, maxs, y = get_min_max_med(data_label)
            mins, maxs, y = np.array(mins), np.array(maxs), np.array(y)
            x = df_parsed["aa_len"][0][:len(y)]
            yerr = np.array([y - mins, maxs - y])
            plt.errorbar(
                x,
                y,
                yerr=yerr,
                label=plot_label,
                # marker=marker,
            )
        
    
    # Plot the time data
    make_plot(time_data_labels)
    plt.title("Execution Time vs. Protein Length")
    plt.xlabel("Protein Length (in amino acids)")
    plt.ylabel("Execution Time (in minutes)")
    plt.legend()
    plt.grid(True)
    plt.yscale(plot_scale)
    plt.savefig("../data/bench_time.pdf", dpi=300)

    # Plot the memory
    make_plot(mem_data_labels)
    plt.title("Memory Usage vs. Protein Length")
    plt.xlabel("Protein Length (in amino acids)")
    plt.ylabel("Memory Usage (in GB)")
    plt.legend()
    plt.grid(True)
    plt.yscale(plot_scale)
    plt.savefig("../data/bench_mem.pdf", dpi=300)


if __name__ == "__main__":
    main()
