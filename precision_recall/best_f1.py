#!/usr/bin/env python3
import pandas as pd
import argparse


import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", default="50mers_F1Scores.csv", type=str, help="name of input csv file")
parser.add_argument("-o", "--output", default="50mers_best_F1_history.csv", type=str, help="name of output csv file")
args = parser.parse_args()

df = pd.read_csv(args.input, sep="\t")

# Find the best F1 for each edit distance
best = df.loc[df.groupby("Edit dist")["F1 Score"].idxmax()]

# Save to a new tab-delimited file
best.to_csv(args.output, sep="\t", index=False,float_format="%.4f")

