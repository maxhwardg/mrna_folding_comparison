import pandas as pd
import matplotlib.pyplot as plt

file = open("../data/bench_runs/tmp", "r")
raw_data = file.read()
plot_scale = "linear"

# Function to parse the raw data
def parse_data(raw_data):
    data = {
        "aa_len": [],
        "lineardesign_time": [],
        "cdsfold_time": [],
        "derna_time": [],
        "mrnafold_time": []
    }
    for line in raw_data.split("\n"):
        if line.startswith("aa_len:"):
            data["aa_len"].append(int(line.split(":")[1].strip()))
        elif line.startswith("lineardesign time(s):"):
            data["lineardesign_time"].append(float(line.split(":")[1].strip()))
        elif line.startswith("cdsfold time(s):"):
            data["cdsfold_time"].append(float(line.split(":")[1].strip()))
        elif line.startswith("derna time(s):"):
            data["derna_time"].append(float(line.split(":")[1].strip()))
        elif line.startswith("mrnafold time(s):"):
            data["mrnafold_time"].append(float(line.split(":")[1].strip()))
            
    def pad_list(lst, length):
        lst += [None] * (length - len(lst))
    
    pad_list(data["lineardesign_time"], len(data["aa_len"]))
    pad_list(data["cdsfold_time"], len(data["aa_len"]))
    pad_list(data["derna_time"], len(data["aa_len"]))
    pad_list(data["mrnafold_time"], len(data["aa_len"]))
    
    
    return data

# Parse the data
parsed_data = parse_data(raw_data)
df_parsed = pd.DataFrame(parsed_data)

print(df_parsed)

# Plotting
plt.figure(figsize=(12, 8))
plt.plot(df_parsed["aa_len"], df_parsed["lineardesign_time"], label="LinearDesign", marker="o")
plt.plot(df_parsed["aa_len"], df_parsed["cdsfold_time"], label="CDSfold", marker="s")
plt.plot(df_parsed["aa_len"], df_parsed["derna_time"], label="deRNA", marker="^")
# plt.plot(df_parsed["aa_len"], df_parsed["mrnafold_time"], label="mRNAfold", marker="d")

plt.title("Execution Time vs. Amino Acid Length")
plt.xlabel("Amino Acid Length (aa_len)")
plt.ylabel("Execution Time (s)")
plt.legend()
plt.grid(True)
plt.yscale(plot_scale)
plt.show()
plt.savefig("../data/bench.png")
