import pandas as pd
import numpy as np
from matplotlib import pyplot as plt



comment = {}
for j in snakemake.input: 
    with open(j) as f:
        for l in f:
            if l.startswith('#batch'):
                i = l[6:].rstrip()
                comment[j+'_cindex'] = []
                comment[j+'_alpha'] = []
                comment[j+'_features'] = []
            if l.startswith('#cindex'):
                comment[j+'_cindex'].append(float(l.split(',')[0][8:]))
                comment[j+'_alpha'].append(float(l.split(',')[1][7:-2]))
            if l.startswith('#fea'):
                comment[j+'_features'].append(l[10:-1])
comment = pd.DataFrame(comment)

cindex = []
for j in snakemake.input:
    cindex += list(comment[j+'_cindex'])
cindex_ = np.quantile(cindex, 0.5)

batches = pd.read_csv(snakemake.input[0], index_col=[0], comment='#')
batches.columns = [snakemake.input[0]]
for i in snakemake.input[1:]:
    f = pd.read_csv(i, index_col=[0], comment='#')
    f.columns = [i]
    batches = pd.concat([batches, f], axis=1)

for j in snakemake.input:
    feature_drop = ','.join(comment[comment[j+'_cindex'] < cindex_][j+'_features']).split(',')
    if len(feature_drop) > 1:
        batches.loc[[x for x in feature_drop if x != ''], j] = 0

# genes: top 80%, transcripts: top 70%
batches[batches==0]=np.nan
batches['MEAN'] = batches.apply(np.mean, axis=1)
batches['Count'] = batches.apply(lambda x: x.notna().sum(), axis=1)-1
features = batches[batches['Count']>np.ceil(max(batches['Count'])*0.1)]
#features = batches[batches['Count']>0]


pd.DataFrame(features.index, columns=['Feature']).to_csv(snakemake.output[0], index=False)
batches.to_csv(snakemake.output[1], index=True)
with open(snakemake.output[2], 'w') as f:
    f.write('#mean:{a},median:{b}\n'.format(a=np.mean(cindex), b=cindex_))
comment.to_csv(snakemake.output[2], index=False, mode='a')













