import pandas as pd
import numpy as np


t = pd.read_csv('transcripts1.csv')
for i in range(2, 5):
    t = pd.merge(t, pd.read_csv('transcripts'+str(i)+'.csv'), on='ID')
    print(i)

t = t.set_index('ID')
t = np.log2(t + 1)
t['MEAN'] = t.apply(np.mean, axis=1)
t['VAR'] = t.apply(np.var, axis=1)
t = t[t['MEAN']>0.01]
t = t.sort_values(by='VAR').iloc[-int(t.shape[0]*0.8):]
t = t.reset_index()
t.to_csv('transcripts_001_08.csv', index=False)
