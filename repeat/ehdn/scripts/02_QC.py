#!/usr/bin/env python3

import argparse
import os
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(description="Summarize and plot norm_num_anc_irrs values.")
    parser.add_argument("--input_dir", default="results/01/", required=True, help="Input directory containing .tsv files")
    parser.add_argument("--output_table", default="results/02/summary.tsv", required=True, help="Path to output summary table (TSV)")
    parser.add_argument("--output_dir", default="results/02/", required=True, help="Directory to save plots")
    return parser.parse_args()


def main():
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    thresholds = [1, 2, 5, 10]
    summary = []

    for file_path in glob(os.path.join(args.input_dir, "*.tsv")):
        df = pd.read_csv(file_path, sep="\t")
        sample = os.path.basename(file_path).replace(".tsv", "")
        row = {"Sample": sample}
        for t in thresholds:
            row[f">={t}"] = (df['norm_num_anc_irrs'] >= t).sum()
        summary.append(row)

    summary_df = pd.DataFrame(summary)
    summary_df.set_index("Sample", inplace=True)
    summary_df.to_csv(args.output_table, sep="\t")

    for t in thresholds:
        col = f">={t}"
        plt.figure(figsize=(10, 6))
        sns.boxplot(y=summary_df[col])
        sns.stripplot(y=summary_df[col], color='red', jitter=True)

        Q1 = summary_df[col].quantile(0.25)
        Q3 = summary_df[col].quantile(0.75)
        IQR = Q3 - Q1
        outliers = summary_df[summary_df[col] > Q3 + 1.5 * IQR]
        for sample, val in outliers[col].items():
            plt.text(0, val, sample, ha='left', va='bottom', fontsize=8)

        plt.title(f"Outliers for norm_num_anc_irrs >= {t}")
        plt.ylabel("Number of rows")
        plot_path = os.path.join(args.output_dir, f"outlier_plot_ge{t}.png")
        plt.savefig(plot_path)
        plt.close()


if __name__ == "__main__":
    main()
