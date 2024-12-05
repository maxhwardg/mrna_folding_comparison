import pandas as pd
import matplotlib.pyplot as plt

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
        "mrnafold_mem": []
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
            
    def pad_list(lst, length):
        lst += [None] * (length - len(lst))
    
    pad_list(data["lineardesign_time"], len(data["aa_len"]))
    pad_list(data["cdsfold_time"], len(data["aa_len"]))
    pad_list(data["derna_time"], len(data["aa_len"]))
    pad_list(data["mrnafold_time"], len(data["aa_len"]))
    pad_list(data["lineardesign_mem"], len(data["aa_len"]))
    pad_list(data["cdsfold_mem"], len(data["aa_len"]))
    pad_list(data["derna_mem"], len(data["aa_len"]))
    pad_list(data["mrnafold_mem"], len(data["aa_len"]))
    
    def convert_to_gb(lst):
        return [None if x is None else x / 1024**3 for x in lst]
    
    data["lineardesign_mem"] = convert_to_gb(data["lineardesign_mem"])
    data["cdsfold_mem"] = convert_to_gb(data["cdsfold_mem"])
    data["derna_mem"] = convert_to_gb(data["derna_mem"])
    data["mrnafold_mem"] = convert_to_gb(data["mrnafold_mem"])
    
    
    return data

# Parse the data
parsed_data = parse_data(raw_data)
df_parsed = pd.DataFrame(parsed_data)

print(df_parsed)

# Plot the time

# Plotting
plt.figure(figsize=(12, 8))
plt.plot(df_parsed["aa_len"], df_parsed["lineardesign_time"], label="LinearDesign", marker="o")
plt.plot(df_parsed["aa_len"], df_parsed["cdsfold_time"], label="CDSfold", marker="s")
plt.plot(df_parsed["aa_len"], df_parsed["derna_time"], label="DERNA", marker="^")
# plt.plot(df_parsed["aa_len"], df_parsed["mrnafold_time"], label="mRNAfold", marker="d")

plt.title("Execution Time vs. Protein Length")
plt.xlabel("Protein Length (in amino acids)")
plt.ylabel("Execution Time (s)")
plt.legend()
plt.grid(True)
plt.yscale(plot_scale)
# plt.show()
plt.savefig("../data/bench_time.pdf", dpi=300)

# Plot the memory
plt.figure(figsize=(12, 8))
plt.plot(df_parsed["aa_len"], df_parsed["lineardesign_mem"], label="LinearDesign", marker="o")
plt.plot(df_parsed["aa_len"], df_parsed["cdsfold_mem"], label="CDSfold", marker="s")
plt.plot(df_parsed["aa_len"], df_parsed["derna_mem"], label="DERNA", marker="^")
# plt.plot(df_parsed["aa_len"], df_parsed["mrnafold_mem"], label="mRNAfold", marker="d")

plt.title("Memory Usage vs. Protein Length")
plt.xlabel("Protein Length (in amino acids)")
plt.ylabel("Memory Usage (in Gigabytes)")
plt.legend()
plt.grid(True)
plt.yscale(plot_scale)
# plt.show()
plt.savefig("../data/bench_mem.pdf", dpi=300)
