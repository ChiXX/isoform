import pandas as pd
import numpy as np
import warnings
from sksurv.linear_model import CoxnetSurvivalAnalysis, CoxPHSurvivalAnalysis
from sksurv.preprocessing import OneHotEncoder
from sksurv.metrics import cumulative_dynamic_auc
from sklearn.model_selection import GridSearchCV, KFold, train_test_split, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from matplotlib import pyplot as plt

def load_data(dataP, spInfo, features):
    data = pd.read_csv(dataP).set_index('ID')
    sample_info = pd.read_csv(spInfo)
#    data = np.log2(data + 1)
#    data['MEAN'] = data.apply(np.mean, axis=1)
#    data['VAR'] = data.apply(np.var, axis=1)
#    data = data[data['MEAN']>0.01]
#    data = data.sort_values(by='VAR').iloc[-int(data.shape[0]*0.8):]
    data.drop(['MEAN','VAR'], axis=1, inplace=True)
    data = data.T.reset_index().rename(columns={'index':'SAMPLE'})
    data.columns.names = [None]
    if isinstance(features, pd.DataFrame):
        data = data[['SAMPLE']+list(features.Feature)]
    if 'RFi' in event:
        data = pd.merge(sample_info[['SAMPLE', 'RFi_days', 'RFi_event', 'RFS_all_event', 'age', 'treatGroup']], data, on='SAMPLE')
        data = data[data.RFi_days.notna()].reset_index(drop=True)
        data = data[data.apply(lambda x: False if (x['RFi_event'] == 0 and x['RFi_days'] < 600) or (x['RFi_event'] == False and x['RFS_all_event'] == True) else True, axis = 1)].reset_index(drop=True)
        data.drop(['RFS_all_event'], axis=1, inplace=True)
    else:
        data = pd.merge(sample_info[['SAMPLE', days, event, 'age', 'treatGroup']], data, on='SAMPLE')
        data = data[data[event].notna()].reset_index(drop=True)
        data = data[data[days].notna()].reset_index(drop=True)
    x, y = data.iloc[::, 3:], data.iloc[::, 1:3]
    treatment = x.treatGroup.astype('category')
    x.drop('treatGroup', axis=1, inplace=True)
    x = pd.DataFrame(StandardScaler().fit_transform(x), columns=x.columns).dropna(axis=1) 
    x['treatGroup'] = treatment
    x = OneHotEncoder().fit_transform(x).replace(np.nan, 0)
    return x, y

def gcv_plots(grids, figPath):
    cv_results = pd.DataFrame(grids.cv_results_)
    alphas = cv_results.param_coxnetsurvivalanalysis__alphas.map(lambda x: x[0])
    mean = cv_results.mean_test_score
    std = cv_results.std_test_score
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(alphas, mean)
    ax.fill_between(alphas, mean - std, mean + std, alpha=.15)
    ax.set_xscale("log")
    ax.set_ylabel("concordance index")
    ax.set_xlabel("alpha")
    ax.axvline(grids.best_params_["coxnetsurvivalanalysis__alphas"], c="C1")
    ax.axhline(0.5, color="grey", linestyle="--")
    ax.grid(True)
    fig.savefig(figPath)
    plt.close()
def cox_linear(x_train, batchid, indexes, threads):
    os_type = {'names':('events', 'time'), 'formats':('?', '<f8')}
    kept = ['age', 'treatGroup=EndoImmu', 'treatGroup=Endo', 'treatGroup=EndoCyto', 'treatGroup=EndoCytoImmu', 'treatGroup=CytoImmu']
    y_train2 = np.array([(m, n) for m, n in zip(ytrain[event].astype(bool), ytrain[days])], dtype=os_type)
    base = CoxPHSurvivalAnalysis()
    base_score = cross_val_score(base, xtrain[kept], y_train2, cv=indexes, n_jobs=threads)
    gcv = GridSearchCV(
        make_pipeline(CoxPHSurvivalAnalysis()),
        param_grid={"coxphsurvivalanalysis__alpha": [0]},
        cv=indexes,
        error_score=np.mean(base_score),
        n_jobs=threads).fit(xtrain[x_train.columns], y_train2)
    best_model = gcv.best_estimator_.named_steps["coxphsurvivalanalysis"]
    selection_coefs = pd.DataFrame(
        best_model.coef_,
        index=x_train.columns,
        columns=[batchid]
    )
    return selection_coefs, np.mean(base_score), gcv 
    
def coefs_batch(x_train, batchid, indexes, threads=8):

    os_type = {'names':('events', 'time'), 'formats':('?', '<f8')}
    #x_train = x_train.iloc[np.append(indexes[0][0], indexes[0][1]),:]
    #y_train = ytrain.iloc[np.append(indexes[0][0], indexes[0][1]),:]
    #y_train = np.array([(m, n) for m, n in zip(y_train[event].astype(bool), y_train[days])], dtype=os_type)
    y_train2 = np.array([(m, n) for m, n in zip(ytrain[event].astype(bool), ytrain[days])], dtype=os_type)

    # test_run = CoxnetSurvivalAnalysis(l1_ratio=0.9, alpha_min_ratio=0.05, n_alphas=30)
    # test_run.fit(x_train, y_train)
    # alphas = 10. ** np.linspace(-2.3, -1.8, 8)
    gcv = GridSearchCV(
        make_pipeline(CoxnetSurvivalAnalysis(l1_ratio=0.9, tol=1e-10, max_iter=10000)),
        # param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in alphas]},
        param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in [0.02, 0.03, 0.04, 0.05, 0.06, 0.07]]},
        cv=indexes,
        error_score=0.6,
        n_jobs=threads).fit(xtrain[x_train.columns], y_train2)
    best_model = gcv.best_estimator_.named_steps["coxnetsurvivalanalysis"]
    best_coefs = pd.DataFrame(
        best_model.coef_,
        index=x_train.columns,
        columns=[batchid]
    )
    return best_coefs.iloc[6:,:], gcv


def training_batch(x_train, y_train, iteration):
    if 'RFi' in event:
        kept = ['age', 'treatGroup=Endo', 'treatGroup=EndoCyto', 'treatGroup=EndoCytoImmu']
    else:
        kept = ['age', 'treatGroup=EndoImmu', 'treatGroup=Endo', 'treatGroup=EndoCyto', 'treatGroup=EndoCytoImmu', 'treatGroup=CytoImmu'] 
    xa = x_train[kept]
    xb = x_train.drop(kept, axis=1)
    yTrue = y_train[y_train[event] == True].reset_index()
    lq = int(yTrue.describe().iloc[4,0]) # lower quartile
    yFalse = y_train[y_train.apply(lambda x: True if x[event] == False and x[days] > lq else False, axis=1)].reset_index()
    if yTrue.shape[0] > yFalse.shape[0]:
        majority = True
    else:
        majority = False
    batches = pd.DataFrame(index=xb.columns)
    
    for i in range(iteration):
        truePart = []
        falsePart = []
        indexes = []
        yTrueKF = KFold(n_splits=5, shuffle=True, random_state=1723).split(yTrue)
        yFalseKF = KFold(n_splits=5, shuffle=True, random_state=1938).split(yFalse)
        for a, j in yTrueKF:
            np.random.shuffle(j) # random part
            if majority:
                truePart.append(yTrue.iloc[j]['index'][:int(yFalse.shape[0]/5)])
            else:
                truePart.append(yTrue.iloc[j]['index'])

        for a, j in yFalseKF:
            np.random.shuffle(j) # random part
            if not majority:
                falsePart.append(yFalse.iloc[j]['index'][:int(yTrue.shape[0]/5)])
            else:
                falsePart.append(yFalse.iloc[j]['index'])
        for a in range(5):
            to_del = [0,1,2,3,4]
            to_del.pop(a)
            trainPart = np.array([], dtype=int)
            valPart = np.append(truePart[a], falsePart[a])
            for j in to_del:
                tep = np.append(truePart[j], falsePart[j])
                trainPart = np.append(trainPart, tep)
            indexes.append([trainPart, valPart])
        if x_train.shape[1]//int(yTrue.shape[0]*1.6) == 0:
            kf = [[np.nan, 0]]
        else:
            kf = KFold(n_splits=x_train.shape[1]//int(yTrue.shape[0]*1.6)+1, shuffle=True).split(xb.T) # splits for features
        batch = pd.DataFrame(columns=['batch'+str(i+1)])
        with open(out, 'a') as f:
            f.write('#batch'+str(i+1)+'\n')
            for train_index, test_index in kf:
                if isinstance(test_index, int):
                    X_sub = xb.T
                else:
                    X_sub = xb.T.iloc[test_index,:]
                xc = pd.concat([xa, X_sub.T], axis=1)
                # xc = xc.iloc[np.append(indexes[0][0], indexes[0][1]),:]
                b, gcv = coefs_batch(xc, 'batch'+str(i+1), indexes, threads=16)
                f.write('#cindex:'+str(gcv.best_score_)+',alpha:'+str(gcv.best_params_["coxnetsurvivalanalysis__alphas"])+'\n')
                # f.write('#baseindex={m}, cindex={n}\n'.format(m=baseScore, n=gcv.best_score_))
                f.write('#features:'+','.join(b[b['batch'+str(i+1)]!=0].index)+'\n')
                # f.write('# '+str(indexes)+'\n')
                batch = pd.concat([batch, b]) # vertical
                # gcv_plots(gcv, '/home/chixu/isoform/denovo/transcripts/figs/'+out.split('/')[-1].split('.')[0]+'_'+str(len(batch.index))+'.jpg')
                #pd.DataFrame(gcv.cv_results_).to_csv('/home/chixu/isoform/normal/genes/figs/'+out.split('/')[-1].split('.')[0]+'_'+str(len(batch.index))+'.csv')
        batches = pd.concat([batches, batch], axis=1) # horizontal
    return batches



event = snakemake.params[0]
days = snakemake.params[1]


if snakemake.input[1] == snakemake.input[2]:
    features = ''
else:
    features = pd.read_csv(snakemake.input[2])    
x, y = load_data(snakemake.input[0], snakemake.input[1], features)
xtrain, xval, ytrain, yval = train_test_split(x, y, test_size=0.2, random_state=954)
out = snakemake.output[0]
xtrain = xtrain.reset_index(drop=True)
ytrain = ytrain.reset_index(drop=True)
batch = training_batch(xtrain, ytrain, iteration=1)
# batch = training_batch(xtrain, ytrain, [21])
batch.to_csv(out, mode='a')

