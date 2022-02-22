import pandas as pd
import numpy as np
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.preprocessing import OneHotEncoder
from sksurv.metrics import cumulative_dynamic_auc
from sklearn.model_selection import GridSearchCV, KFold, train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler


event = snakemake.params[0]
days = snakemake.params[1]


def load_data(dataP, spInfo):
    data = pd.read_csv(dataP).set_index('ID')
    sample_info = pd.read_csv(spInfo)
    data.drop(['MEAN','VAR'], axis=1, inplace=True)
    data = data.T.reset_index().rename(columns={'index':'SAMPLE'})
    data.columns.names = [None]
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

x, y = load_data(snakemake.input[0], snakemake.input[1])
os_type = {'names':('events', 'time'), 'formats':('?', '<f8')}
y = np.array([(a, b) for a, b in zip(y[event].astype(bool), y[days])], dtype=os_type)

xtrain, xval, ytrain, yval = train_test_split(x, y, test_size=0.2, random_state=954)

templete = CoxnetSurvivalAnalysis(l1_ratio=0.9, alpha_min_ratio=0.01, n_alphas=20)
templete.fit(xtrain, ytrain)
cv = KFold(n_splits=5, shuffle=True, random_state=0)
grids = GridSearchCV(
    make_pipeline(CoxnetSurvivalAnalysis(l1_ratio=0.9)),
    param_grid={"coxnetsurvivalanalysis__alphas": [[v] for v in templete.alphas_]},
    cv=cv,
    error_score=0.5,
    n_jobs=16).fit(xtrain, ytrain)

best_model = grids.best_estimator_.named_steps["coxnetsurvivalanalysis"]
best_coefs = pd.DataFrame(
    best_model.coef_,
    index=xtrain.columns,
    columns=['coefficient']
)

pd.DataFrame(grids.cv_results_).to_csv(snakemake.output[0], index=False)
best_coefs.query("coefficient != 0").to_csv(snakemake.output[1])
