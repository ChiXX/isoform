import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split


uni = pd.read_csv(snakemake.input[0], comment='#')
onerun = pd.read_csv(snakemake.input[1], comment='#')
ens = pd.read_csv(snakemake.input[2], comment='#')
sample_info = pd.read_csv(snakemake.input[3])

fo = snakemake.output[0] 

event = snakemake.params[0]
days = snakemake.params[1]

with open(snakemake.input[0]) as fi:
    for l in fi:
        if l.startswith('# base_mean_auc:'):
            base_mean = float(l.split(':')[1].strip())

label = os.path.basename(fo).split('-')

if 'gs' in label[1]:
    data_set = 'gene'
else:
    data_set = 'transcript'

if 'normal' in label[0]:
    version = 'V32'
elif 'denovo' in label[0]:
    version = 'de-novo'
else:
    version = 'V27'

if 'ki67_censored' in label[1]:
    label = 'Ki67_censored'
elif 'ki67_uncensored' in label[1]:
    label = 'Ki67_uncensored'
else:
    label = event.split("_")[0]


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



fg = plt.figure(figsize=(15, 5), dpi=128)
ax = fg.add_subplot(111)

ax.plot(uni.time, uni.AUC, marker=".", color='darkorange', label='Univariate analysis')
ax.axhline(uni['mean'][0], linestyle="--", color='bisque', label='Accumulated AUC = %.3f'%uni['mean'][0])

ax.plot(onerun.time, onerun.AUC, marker=".", color='dodgerblue', label='One run analysis')
ax.axhline(onerun['mean'][0], linestyle="--", color='lightskyblue', label='Accumulated AUC = %.3f'%onerun['mean'][0])

ax.plot(ens.time, ens.AUC, marker=".", color='forestgreen', label='Ensemble method')
ax.axhline(ens['mean'][0], linestyle="--", color='palegreen', label='Accumulated AUC = %.3f'%ens['mean'][0])

ax.plot(ens.time, ens.base, marker=".", color='darkgrey', label='Baseline')
ax.axhline(base_mean, linestyle="--", color='silver', label='Accumulated AUC = %.3f'%base_mean)

plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

ax.set_xticklabels([]) 
ax.set_ylabel("Time-Dependent AUC", fontsize=18)
# ax.set_title(f'{label} label using {data_set} set assembled by {version} annotation file', fontsize=22)
plt.legend(fontsize=12)

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
