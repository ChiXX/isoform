import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

path = os.path.split(snakemake.output[0])[0]
label = os.path.split(snakemake.output[0])[1].split('-')[0]

N_round = [i for i in range(1,10)]
dataset = [path for i in range(9)]
performance = []
features = []

for i in N_round:
    with open(os.path.join(path, label+f'-describe.Round{i}')) as f:
        performance.append(float(f.readline().split(',')[0][6:]))
    features.append(pd.read_csv(os.path.join(path, label+f'-features.Round{i}')).shape[0])

pd.read_csv(os.path.join(path, label+f'-features.Round{performance.index(max(performance))+1}')).to_csv(snakemake.output[0], index=False)
pd.DataFrame({'dataset':dataset, 'N_round':N_round, 'performance':performance, 'features':features}).to_csv(snakemake.output[1], index=False)

