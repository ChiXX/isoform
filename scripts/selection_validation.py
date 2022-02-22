import pandas as pd
import numpy as np
import warnings
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.preprocessing import OneHotEncoder
from sksurv.metrics import cumulative_dynamic_auc
from sklearn.model_selection import GridSearchCV, KFold, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler

def load_data(dataP, spInfo, features):
    data = pd.read_csv(dataP).set_index('ID')
    sample_info = pd.read_csv(spInfo)
    sample_info['Ki67_event'] = True
    data.drop(['MEAN','VAR'], axis=1, inplace=True)
    data = data.T.reset_index().rename(columns={'index':'SAMPLE'})
    data.columns.names = [None]
    #data = data[['SAMPLE']+list(features)]
    data = pd.merge(sample_info[['SAMPLE', 'Ki67_fake_days', 'Ki67_event', 'age', 'treatGroup', 'tumor_size']], data, on='SAMPLE')
    data = data[data.Ki67_fake_days.notna()].reset_index(drop=True)
    data = data[data.tumor_size.notna()].reset_index(drop=True)
    x, y = data.iloc[::, 3:], data.iloc[::, 1:3]
    treatment = x.treatGroup.astype('category')
    x.drop('treatGroup', axis=1, inplace=True)
    x = pd.DataFrame(StandardScaler().fit_transform(x), columns=x.columns).dropna(axis=1) 
    x['treatGroup'] = treatment
    x = OneHotEncoder().fit_transform(x).replace(np.nan, 0)
    return x, y

def coefs_batch(x_train, y_train, batch, threads=8):
    test_run = CoxnetSurvivalAnalysis(l1_ratio=0.9, alpha_min_ratio=0.01, n_alphas=50)
    test_run.fit(x_train, y_train)
    if len(test_run.alphas_) == 50:
        test_alphas = test_run.alphas_[15:-5]
    else:
        test_alphas = test_run.alphas_[:]
    cv = KFold(n_splits=5, shuffle=True, random_state=0)
    gcv = GridSearchCV(
        make_pipeline(CoxnetSurvivalAnalysis(l1_ratio=0.9)),
        param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in test_alphas]},
        cv=cv,
        error_score=0.5,
        n_jobs=threads).fit(x_train, y_train)
    best_model = gcv.best_estimator_.named_steps["coxnetsurvivalanalysis"]
    best_coefs = pd.DataFrame(
        best_model.coef_,
        index=x_train.columns,
        columns=[batch]
    )
    return best_coefs.iloc[7:,:], gcv.best_score_, gcv.best_params_["coxnetsurvivalanalysis__alphas"]

def training_batch(x_train, y_train, seeds):
    kept = ['age', 'tumor_size', 'treatGroup=EndoImmu', 'treatGroup=Endo', 'treatGroup=EndoCyto', 'treatGroup=EndoCytoImmu', 'treatGroup=CytoImmu']
    xa = x_train[kept]
    xb = x_train.drop(kept, axis=1)
    os_type = {'names':('events', 'time'), 'formats':('?', '<f8')}
    batches = pd.DataFrame(index=xb.columns)
    i = 1
    for s in seeds:
        #if x_train.shape[1]//int(336*1.6) == 0:
        #    kf = [[np.nan, 0]]
        #else:
        kf = KFold(n_splits=x_train.shape[1]//int(336*1.6)+1, shuffle=True, random_state=s).split(xb.T)
        np.random.seed(s)
        batch = pd.DataFrame(columns=['batch'+str(i)])
        with open(out, 'a') as f:
            f.write('#batch'+str(i)+'\n')
            for train_index, test_index in kf:
                #if test_index == 0:
                #    X_sub = xb.T
                #else:
                X_sub = xb.T.iloc[test_index,:]
                state = np.random.randint(99999)
                xc = pd.concat([xa, X_sub.T], axis=1)
                yTrue = y_train.sample(n=336, random_state=state)
                yFalse = y_train.drop(yTrue.index, axis=0).sort_values(by='Ki67_fake_days', ascending=False).iloc[:1000,:].sample(n=336, random_state=state)
                yFalse['Ki67_fake_days'] = yFalse['Ki67_fake_days']*0.1
                yFalse['Ki67_event'] = False
                index = yTrue.index.append(yFalse.index)
                y_train2 = yTrue.append(yFalse).sample(frac=1, random_state=state)
                y_train2 = np.array([(m, n) for m, n in zip(y_train2['Ki67_event'].astype(bool), y_train2['Ki67_fake_days'])], dtype=os_type)
                x_train2 = xc.loc[index,:].sample(frac=1, random_state=state)
                b, score, alpha = coefs_batch(x_train2, y_train2, 'batch'+str(i), threads=16)
                f.write('#cindex:'+str(score)+',alpha:'+str(alpha)+'\n')
                f.write('#features:'+','.join(b[b['batch'+str(i)]!=0].index)+'\n')
                batch = pd.concat([batch, b]) # vertical
        i += 1
        batches = pd.concat([batches, batch], axis=1) # horizontal
    return batches


x, y = load_data(snakemake.input[0], snakemake.input[1], [])
#features = pd.read_csv(snakemake.input[2])    
#x, y = load_data(snakemake.input[0], snakemake.input[1], features.Feature)
xtrain, xval, ytrain, yval = train_test_split(x, y, test_size=0.2, random_state=954)
out = snakemake.output[0]
np.random.seed(int(out.split('/')[-1].split('.')[0]))
batch = training_batch(xtrain, ytrain, np.random.randint(99999, size=12))
# batch = training_batch(xtrain, ytrain, [21])
batch.to_csv(out, mode='a')

