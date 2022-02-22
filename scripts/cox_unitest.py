from scipy import stats
from lifelines import CoxPHFitter
import pandas as pd
import numpy as np

event = snakemake.params[0]
days = snakemake.params[1]

def survival(df, t):
    try:
        return float(CoxPHFitter().fit(df[[t, days, event]], days, event_col=event).summary['p'])
    except:
        return np.nan

data = pd.read_csv(snakemake.input[0]).set_index('ID')
data = data.T.reset_index().rename(columns={'index':'SAMPLE'})
data.columns.names = [None]

if 'RFi' in event:
    data = pd.merge(pd.read_csv(snakemake.input[1])[['SAMPLE', 'RFi_days', 'RFi_event', 'RFS_all_event']], data, on='SAMPLE')
    data = data[data.RFi_days.notna()].reset_index(drop=True)
    data = data[data.apply(lambda x: False if (x['RFi_event'] == 0 and x['RFi_days'] < 600) or (x['RFi_event'] == False and x['RFS_all_event'] == True) else True, axis = 1)].reset_index(drop=True)
    data.drop(['RFS_all_event'], axis=1, inplace=True)
else:
    data = pd.merge(pd.read_csv(snakemake.input[1])[['SAMPLE', days, event]], data, on='SAMPLE')
    data = data[data[event].notna()].reset_index(drop=True)
    data = data[data[days].notna()].reset_index(drop=True)

ts = data.columns[3:]
ts_sig = {}
for t in ts:
    ts_sig[t] = survival(data, t)

ts_df = pd.DataFrame.from_dict(ts_sig, orient='index', columns=['p-val']).reset_index().rename(columns={'index':'fetures'})
ts_df.to_csv(snakemake.output[0], index=False)
