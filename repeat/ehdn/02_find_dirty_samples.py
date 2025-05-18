#!/usr/bin/env python3


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from glob import glob


thresholds = [1, 2, 5, 10]
summary = []

for file_path in glob("result_step1/*.tsv"):  # ファイルパスを調整
    df = pd.read_csv(file_path, sep="\t")
    sample = os.path.basename(file_path).replace(".tsv", "")
    row = {"Sample": sample}
    for t in thresholds:
        row[f">={t}"] = (df['norm_num_anc_irrs'] >= t).sum()
    summary.append(row)

summary_df = pd.DataFrame(summary)
summary_df.set_index("Sample", inplace=True)
summary_df.to_csv("summary_table.tsv", sep="\t")


for t in thresholds:
    col = f">={t}"
    plt.figure(figsize=(10, 6))
    sns.boxplot(y=summary_df[col])
    sns.stripplot(y=summary_df[col], color='red', jitter=True)

    # 外れ値を注釈
    Q1 = summary_df[col].quantile(0.25)
    Q3 = summary_df[col].quantile(0.75)
    IQR = Q3 - Q1
    outliers = summary_df[summary_df[col] > Q3 + 1.5 * IQR]
    for sample, val in outliers[col].items():
        plt.text(0, val, sample, ha='left', va='bottom', fontsize=8)

    plt.title(f"Outliers for norm_num_anc_irrs >= {t}")
    plt.ylabel("Number of rows")
    plt.savefig(f"outlier_plot_ge{t}.png")
    plt.close()

