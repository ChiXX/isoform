import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split

sample_info = pd.read_csv(snakemake.input[-1])

fo = snakemake.output[0]

event = snakemake.params[0]
days = snakemake.params[1]

label = event.split("_")[0].lower()

de_gs_uni = pd.read_csv(f'denovo/{label}_age_trm_gs/cox_{event.split("_")[0]}_uni_analysis.txt', comment='#')
de_gs_1run = pd.read_csv(f'denovo/{label}_age_trm_gs/cox_{event.split("_")[0]}_1run_analysis.txt', comment='#')
de_gs_ens = pd.read_csv(f'denovo/{label}_age_trm_gs/cox_{event.split("_")[0]}_ens_analysis.txt', comment='#')

v32_gs_uni = pd.read_csv(f'normal/{label}_age_trm_gs/cox_{event.split("_")[0]}_uni_analysis.txt', comment='#')
v32_gs_1run = pd.read_csv(f'normal/{label}_age_trm_gs/cox_{event.split("_")[0]}_1run_analysis.txt', comment='#')
v32_gs_ens = pd.read_csv(f'normal/{label}_age_trm_gs/cox_{event.split("_")[0]}_ens_analysis.txt', comment='#')

v27_gs_uni = pd.read_csv(f'v27/{label}_age_trm_gs/cox_{event.split("_")[0]}_uni_analysis.txt', comment='#')
v27_gs_1run = pd.read_csv(f'v27/{label}_age_trm_gs/cox_{event.split("_")[0]}_1run_analysis.txt', comment='#')
v27_gs_ens = pd.read_csv(f'v27/{label}_age_trm_gs/cox_{event.split("_")[0]}_ens_analysis.txt', comment='#')

de_ts_uni = pd.read_csv(f'denovo/{label}_age_trm_ts/cox_{event.split("_")[0]}_uni_analysis.txt', comment='#')
de_ts_1run = pd.read_csv(f'denovo/{label}_age_trm_ts/cox_{event.split("_")[0]}_1run_analysis.txt', comment='#')
de_ts_ens = pd.read_csv(f'denovo/{label}_age_trm_ts/cox_{event.split("_")[0]}_ens_analysis.txt', comment='#')

v32_ts_uni = pd.read_csv(f'normal/{label}_age_trm_ts/cox_{event.split("_")[0]}_uni_analysis.txt', comment='#')
v32_ts_1run = pd.read_csv(f'normal/{label}_age_trm_ts/cox_{event.split("_")[0]}_1run_analysis.txt', comment='#')
v32_ts_ens = pd.read_csv(f'normal/{label}_age_trm_ts/cox_{event.split("_")[0]}_ens_analysis.txt', comment='#')


if 'RFi' in event:
    sample_info = sample_info[['RFi_days', 'RFi_event', 'RFS_all_event']]
    sample_info = sample_info[sample_info.RFi_days.notna()].reset_index(drop=True)
    sample_info = sample_info[sample_info.apply(lambda x: False if (x['RFi_event'] == 0 and x['RFi_days'] < 600) or (x['RFi_event'] == False and x['RFS_all_event'] == True) else True, axis = 1)].reset_index(drop=True)
    sample_info.drop(['RFS_all_event'], axis=1, inplace=True)
else:
    sample_info = sample_info[[days, event]]
    sample_info = sample_info[sample_info[event].notna()].reset_index(drop=True)
    sample_info = sample_info[sample_info[days].notna()].reset_index(drop=True)

_, test = train_test_split(sample_info, test_size=0.2, random_state=954)
va_times = np.arange(np.percentile(sample_info[days], [5, 95])[0], np.percentile(sample_info[days], [5, 95])[1], 60)

os_type = {'names':('events', 'time'), 'formats':('?', '<f8')}
test = np.array([(a, b) for a, b in zip(test[event].astype(bool), test[days])], dtype=os_type)

fg = plt.figure(figsize=(18, 7), dpi=128)
ax = fg.add_subplot(111)

ax.plot(de_gs_uni.time, de_gs_uni.AUC, marker=".", color='lightcoral',label='de_gs_uni')
ax.plot(de_gs_1run.time, de_gs_1run.AUC, marker=".", color='firebrick',label='de_gs_1run')
ax.plot(de_gs_ens.time, de_gs_ens.AUC, marker=".", color='darkred', label='de_gs_ens')

ax.plot(v32_gs_uni.time, v32_gs_uni.AUC, marker=".", color='chocolate', label='v32_gs_uni')
ax.plot(v32_gs_1run.time, v32_gs_1run.AUC, marker=".", color='peru', label='v32_gs_1run')
ax.plot(v32_gs_ens.time, v32_gs_ens.AUC, marker=".", color='darkorange',label='v32_gs_ens')

ax.plot(v27_gs_uni.time, v27_gs_uni.AUC, marker=".", color='darkseagreen', label='v27_gs_uni')
ax.plot(v27_gs_1run.time, v27_gs_1run.AUC, marker=".", color='forestgreen', label='v27_gs_1run')
ax.plot(v27_gs_ens.time, v27_gs_ens.AUC, marker=".", color='seagreen', label='v27_gs_ens')

ax.plot(de_ts_uni.time, de_ts_uni.AUC, marker=".", color='teal', label='de_ts_uni')
ax.plot(de_ts_1run.time, de_ts_1run.AUC, marker=".", color='skyblue', label='de_ts_1run')
ax.plot(de_ts_ens.time, de_ts_ens.AUC, marker=".", color='steelblue', label='de_ts_ens')

ax.plot(v32_ts_uni.time, v32_ts_uni.AUC, marker=".", color='navy',label='v32_ts_uni')
ax.plot(v32_ts_1run.time, v32_ts_1run.AUC, marker=".", color='slateblue', label='v32_ts_1run')
ax.plot(v32_ts_ens.time, v32_ts_ens.AUC, marker=".", color='purple', label='v32_ts_ens')

ax.plot(v32_ts_ens.time, v32_ts_ens.base, marker=".", color='silver', label='Baseline')

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

ax.set_xticklabels([])
ax.set_ylabel("Time-Dependent AUC", fontsize=18)
plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0, fontsize=12)


cum_events = []
cum_no_events = []
for i_time in va_times:
    i_events = len([i_event for i_event in test if (i_event[0] and i_event[1] <= i_time)])
    i_no_events = len([i_event for i_event in test if (i_event[1] > i_time)])
    cum_events.append(i_events)
    cum_no_events.append(i_no_events)

plt.table(cellText=[cum_events, cum_no_events],
          rowLabels=['Events (Cumulative)', 'No Events (-Cumulative)'],
          colLabels=va_times,
          cellLoc='center',
          loc='bottom')


plt.savefig(fo)
































