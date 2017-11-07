'''Script for identifying QTLs using an ensemble L1 norm'''
import pandas as pd
import numpy as np
from sklearn.preprocessing import Imputer
from sklearn.utils import resample
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LassoCV

def impute_data(X):
    '''
    Prepare input SNPs data by imputing missing data.
    This approach will impute data using the most frequent call for each marker
    '''
    X_c = np.copy(X)
    imputer = Imputer(strategy='most_frequent', axis=0)
    X_imp = imputer.fit_transform(X_c)
    return X_imp

def optimize_weights(X, y, cv):
    # for now just optimize for regression
    # for Lasso regression you need to specify the alphas
    C_values = np.logspace(-5, 4, 10)
    alphas =1/C_values
    model = LassoCV(alphas=alphas, random_state=42)
    model.fit(X,y)
    return model.coef_

def mlQTLs_analysis(X, y, model, n_estimators=1000,
                    max_samples = 0.8, cv = 10, random_state=42):
    '''main function that will call the other functions'''
    all_scores = []
    x=random_state
    for i in range(n_estimators):
        # loop over the number of estimators
        X_res, y_res = resample(X, y, n_samples=np.ceil(X.shape[0]*max_samples),
                                replace=True, random_state=x+1)
        w = optimize_weights(X_res, y_res, cv=cv)
        all_w.append(w)
    # now you will need to:
        #1) Calculate the stability score for each feature
        #2) Sum the scores
    all_w = np.array(all_w)
    sum_w = np.sum(all_w, axis=0)
    bin_matrix = np.where(all_w>0, 1, 0)
    s_scores = np.sum(bin_matrix, axis=1)/n_estimators
