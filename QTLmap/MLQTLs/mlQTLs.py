'''Script for identifying QTLs using an ensemble L1 norm'''
import pandas as pd
import numpy as np
from sklearn.preprocessing import Imputer
from sklearn.utils import resample
from sklearn.model_selection import GridSearchCV
from sklearn.svm import LinearSVR

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
    param = {'C':np.logspace(-5, 4, 10)}
    svr = LinearSVR(random_state=42)
    model = GridSearchCV(svr, param_grid = param,
                         scoring='neg_mean_squared_error', cv=10)
    model.fit(X,y)
    model_final = model.best_estimator_
    model_final.fit(X,y)
    return model_final.coef_

def mlQTLs_analysis(X, y, n_estimators=1000,
                    max_samples = 0.8, cv = 10, random_state=42):
    '''main function that will call the other functions'''
    all_w = []
    for i in range(n_estimators):
        # loop over the number of estimators
        X_res, y_res = resample(X, y,
                                n_samples=int(np.ceil(X.shape[0]*max_samples)),
                                replace=True, random_state=i*random_state)
        w = optimize_weights(X_res, y_res, cv=cv)
        all_w.append(w)
    # now you will need to:
        #1) Calculate the stability score for each feature
        #2) Sum the scores
        #3) for optimization purpose I will return the full  coeff_array
    all_w = np.array(all_w)
    sum_w = np.sum(all_w, axis=0)
    bin_matrix = np.where(all_w>0, 1, 0)
    s_scores = np.sum(bin_matrix, axis=0)
    return all_w, s_scores
