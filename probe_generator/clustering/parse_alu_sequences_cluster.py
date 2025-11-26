import os
import argparse
from collections import defaultdict, Counter
from cluster_data import cluster
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm

parser = argparse.ArgumentParser()
parser.add_argument("--method", choices=["bin", "kde"], help="Choose 'bin' for histogram or 'kde' for kernel density estimation (default: bin)")
parser.add_argument("-i", help="input file location with sequences in fasta format (not yet parsed)")
parser.add_argument("-o", default="alu_sequences.fasta", help="name of the output file in which the parsed sequences will be stored")
parser.add_argument("-m", default="alu-metadata.csv", help="name of the file in which the metadata will be stored")

args = parser.parse_args()

input_file = args.i
output_file = args.o
metadata_file = args.m

alu_subfamilies = defaultdict(list)
metadata = defaultdict(list)
alus = []

with open(input_file) as infile, open('alu_sequences_unfiltered.fasta', 'w') as outfile:
    for line in infile:
        if line.startswith(">"):
            header = line.strip()
            # Check if this is an Alu
            if "|Alu|" in header:
                # Find the ME line (Alu sequence)
                for _ in range(3):
                    next(infile, "").strip()

                me_line = next(infile, "").strip().upper()

                # Only adding alus which have a subfamily alu
                alu_subfamily = header.split("|")[5]
                if alu_subfamily.startswith("Alu"):
                    chrm, pos = header.split("|")[3].split(":")
                    start, end = pos.split("-")    
                    length = int(end) - int(start)
                    alus.append((header, me_line, chrm, alu_subfamily, start, end, length))
                    outfile.write(f"{header}\n")
                    outfile.write(f"{me_line}\n")

if not os.path.exists("clusters.txt"):
    # Clusters alu sequences into similar groups
    max_cluster_id = cluster()

# More than one cluster can be used, in this case we choose the cluster with the max amount of sequences
CLUSTERS = [max_cluster_id]
cluster = []
# Holds the alu sequences pertaining to CLUSTERS
for CLUSTER in CLUSTERS:
    with open("clusters.txt") as infile:
        cluster_found = False
        for line in infile:
            if CLUSTER in line:
                cluster_found = True
                continue

            if cluster_found:
                if "Cluster" in line:
                    break
                cluster.append(line.strip())

    # removing empty entry
    cluster = cluster[:-1]

alu_lengths = []
with open(output_file, "w") as outfile:
    for header, me_line, chrm, alu_subfamily, start, end, length in alus:
        if me_line in cluster:
            outfile.write(f"{header}\n{me_line}\n")
            metadata[chrm].append((alu_subfamily, start, end))
            alu_lengths.append(length)

with open(metadata_file, "w") as outfile:
    outfile.write("Chr,Start,End,SubFamily\n")
    for chrm in metadata.keys():
        for name, start, end in metadata[chrm]:
            outfile.write(f"{chrm[3:]},{start},{end},{name}\n")

def plot_freq_dist():
    frequency = Counter(alu_lengths)
    
    # Prepare histogram data
    bin_labels = sorted(frequency.keys())
    frequencies = [frequency[length] for length in bin_labels]
    
    # Create figure
    plt.figure(figsize=(12, 6))
    
    # Plot histogram bars
    x_pos = np.arange(len(bin_labels))
    plt.bar(x_pos, frequencies, width=0.6, alpha=0.7, 
            color='steelblue', label='Observed Frequency')
    plt.xticks(x_pos, bin_labels, rotation=45)
    
    # Calculate normal distribution parameters from raw data
    mu, std = np.mean(alu_lengths), np.std(alu_lengths)
    
    # Create matching x-values for normal curve
    x_normal = np.linspace(min(alu_lengths), max(alu_lengths), 100)
    
    # Calculate normal curve - scale to match histogram counts
    # The scaling factor is: total_counts * bin_width
    bin_width = 1  # Assuming each bar represents 1bp difference
    scale_factor = len(alu_lengths) * bin_width
    normal_curve = norm.pdf(x_normal, mu, std) * scale_factor
    
    # Plot normal curve using the same x-axis positions
    # We need to map the x_normal values to their corresponding bin positions
    x_plot_positions = (x_normal - min(bin_labels)) * (len(bin_labels)-1)/(max(bin_labels)-min(bin_labels))
    
    plt.plot(x_plot_positions, normal_curve, 'r-', linewidth=2,
             label=f'Normal Fit (μ={mu:.1f}, σ={std:.1f})')
    
    # Format plot
    plt.xlabel('Alu Sequence Length (bp)')
    plt.ylabel('Frequency')
    plt.title('Alu Length Distribution with Normal Fit')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

def plot_dist_KDE():
    # Your data
    data = alu_lengths

    # Calculate normal distribution
    mu, std = np.mean(data), np.std(data)
    x = np.linspace(min(data), max(data), 100)
    p = norm.pdf(x, mu, std)

    # Plot KDE (no binning involved)
    plt.figure(figsize=(10, 5))
    sns.kdeplot(data, fill=True, label='KDE (No binning)', linewidth=2)

    # Overlay normal curve
    plt.plot(x, p, 'r--', label='Normal Distribution')

    plt.xlabel('Value')
    plt.ylabel('Density')
    plt.title('KDE with Normal Distribution Overlay')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if args.method == "bin":
    plot_freq_dist()

elif args.method == "kde":
    plot_dist_KDE()
