import pandas as pd
import numpy as np
import os
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.preprocessing import OneHotEncoder
from sksurv.metrics import (
    as_concordance_index_ipcw_scorer,
    as_cumulative_dynamic_auc_scorer,
    as_integrated_brier_score_scorer,
    cumulative_dynamic_auc
)
from sklearn.model_selection import GridSearchCV, KFold, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler


event = snakemake.params[0]
days = snakemake.params[1]


def load_data(dataP, spInfo, features):
    data = pd.read_csv(dataP).set_index('ID')
    sample_info = pd.read_csv(spInfo)
    data.drop(['MEAN','VAR'], axis=1, inplace=True)
    data = data.T.reset_index().rename(columns={'index':'SAMPLE'})
    data.columns.names = [None]
    data = data[['SAMPLE']+features]
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

output = snakemake.output[0]
method = output.split('_')[-2]

if method == 'uni':
    features = pd.read_csv(os.path.join(os.path.dirname(output), f'{event.split("_")[0]}_uni_fdr.csv'))
    features = features[features['p.adjust']<0.05]
    features = features['fetures'].to_list()
elif method == 'ens':
    features = pd.read_csv(os.path.join(os.path.dirname(output), f'{event.split("_")[0]}-best.feature'))
    features = features['Feature'].to_list()
elif method == '1run':
    which = output.split('/')[0]
    dataset = lambda x: 'genes' if 'gs' in x else 'transcripts'
    features = pd.read_csv(os.path.join(os.path.dirname(output), f'{event.split("_")[0]}_1run_{which}_{dataset(output.split("/")[1])}_best_coefs.csv'), index_col=0).reset_index()
    features = features['index'].to_list()
    features = [x for x in features if x != 'age' and not x.startswith('treat')]
else:
    features = []

x, y = load_data(snakemake.input[0], snakemake.input[1], features)
os_type = {'names':('events', 'time'), 'formats':('?', '<f8')}
y = np.array([(a, b) for a, b in zip(y[event].astype(bool), y[days])], dtype=os_type)

xtrain, xval, ytrain, yval = train_test_split(x, y, test_size=0.2, random_state=954)

lower, upper = np.percentile(pd.DataFrame(ytrain).time, [25, 75])
grids_times = np.arange(lower, upper + 1)

cv = KFold(n_splits=3, shuffle=True, random_state=0)
if method == 'uni':
    alphas = 10. ** np.linspace(-1.8, -1, 5)
else:
    alphas = 10. ** np.linspace(-1.9, -1.5, 5)

grids_cindex = GridSearchCV(
    as_concordance_index_ipcw_scorer(CoxnetSurvivalAnalysis(l1_ratio=0.9, tol=1e-15, max_iter=10000), tau=grids_times[-1]),
    param_grid={"estimator__alphas": [[v] for v in alphas]},
    cv=cv,
    error_score=0.5,
    n_jobs=16).fit(xtrain, ytrain)

cph = CoxnetSurvivalAnalysis(l1_ratio=0.9, alphas=[grids_cindex.best_params_["estimator__alphas"]])
cph.fit(xtrain, ytrain)

best_coefs = pd.DataFrame(cph.coef_, index=xtrain.columns, columns=['coefficient'])
non_zero = np.sum(best_coefs.iloc[:, 0] != 0)
non_zero_coefs = best_coefs.query("coefficient != 0")

va_times = np.arange(np.percentile(pd.DataFrame(yval).time, [5, 95])[0], np.percentile(pd.DataFrame(yval).time, [5, 95])[1], 60)
cph_risk_scores = cph.predict(xval)
auc, mean_auc = cumulative_dynamic_auc(ytrain, yval, cph_risk_scores, va_times)



baseline_feature = [x for x in xtrain.columns if x == 'age' or x.startswith('treat')]
base = CoxnetSurvivalAnalysis(l1_ratio=0.9, tol=1e-20, max_iter=10000, alphas=[0])
base.fit(xtrain[baseline_feature], ytrain)
base_risk_scores = base.predict(xval[baseline_feature])
base_auc, base_mean = cumulative_dynamic_auc(ytrain, yval, base_risk_scores, va_times)


with open(output, 'w') as fo:
    fo.write(f'## trainnig \n')
    fo.write(f'# mean cindex: {grids_cindex.cv_results_["mean_test_score"]} \n')
    fo.write(f'# best alpha: {grids_cindex.best_params_["estimator__alphas"]} \n')
    fo.write(f'## validation \n')
    fo.write(f'# cindex:{cph.score(xval, yval)} \n')
    fo.write(f'# features in: {len(features)}+age+treatment, features out: {non_zero} \n')
    fo.write(f'## baseline \n')
    fo.write(f'## features: {baseline_feature} \n')
    fo.write(f'# cindex: {base.score(xval[baseline_feature], yval)} \n')
    fo.write(f'# base_mean_auc: {base_mean} \n')


pd.DataFrame({'AUC':auc, 'base':base_auc, 'time':va_times, 'mean':[mean_auc]*len(auc)}).to_csv(output, index=False, mode='a')



