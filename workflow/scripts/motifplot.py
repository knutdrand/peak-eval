import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

get_name = lambda parts: f"{parts[0]}:{parts[1]}-{parts[2]}"
matches = {line.split("\t")[2] for line in open(snakemake.input.matches) if not line.startswith("#") and line.strip()}
hits = [get_name(line.split()) in matches for line in open(snakemake.input.peaks)]
ratio = np.cumsum(hits)/np.arange(1, len(hits)+1)
plt.plot(ratio)
plt.savefig(snakemake.output[0])
pd.DataFrame({"x": np.arange(ratio.size), "y": ratio}).to_pickle(snakemake.output[1])
