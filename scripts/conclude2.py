import pandas as pd

pd.concat([pd.read_csv(snakemake.input[i]) for i in range(4)]).to_csv(snakemake.output[0], index=False)
