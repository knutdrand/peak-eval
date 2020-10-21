from pathlib import PurePath
import seaborn as sns
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

filenames = snakemake.input
plt.figure(figsize=(10, 10))
dfs = [pd.read_pickle(filename) for filename in filenames]
xmax = min(df["x"].max() for df in dfs)
for df, filename in zip(dfs, filenames):
    df["name"] = filename.split("/")[1]
    print(filename, df)

df = pd.concat(dfs)
p = sns.lineplot(data=df, x="x", y="y", hue="name")
plt.xlabel("Rank")
plt.ylabel("Proportion motif hits")
plt.xlim((0, xmax))
plt.title("Motif_Plot")
plt.savefig(snakemake.output[0])
