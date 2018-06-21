#!/usr/bin/env python2

'''
Input: (1) Path to db.
Output files are created in cwd: (1) ML_scores_output.csv: stores Neg mean squared error values for each predictor/combination.
(2) ML_pred_vs_true_dipeptide.csv: stores predicted OGT and true OGT for dipeptide frequencies.
(3) ML_pred_vs_true_combined_wo_dipep.csv: stores predicted OGT and true OGT for the combination of all features except dipeptide frequencies.
'''

from datetime import datetime
import sqlite3
import numpy as np
from numpy import ravel
from sklearn.model_selection import cross_val_score, cross_val_predict
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from sklearn.svm import SVR
from sklearn.preprocessing import Imputer
from sklearn import preprocessing

db = sys.argv[1] 
conn = sqlite3.connect(db)
c = conn.cursor()
scoreType = 'neg_mean_squared_error'
outputFile = './ML_scores_output.csv'
predFile_dipep = './ML_pred_vs_true_dipeptide.csv'
predFile_combined_wo_dipep = './ML_pred_vs_true_combined_wo_dipep.csv'

def FeaturePrediction(X, y, feature_name, data_Dict):

    imp = Imputer(missing_values='NaN', strategy='mean', axis=0)

    # Impute missing data
    imp.fit(X)
    X_transformed = imp.transform(X)
    scaler = preprocessing.StandardScaler().fit(X_transformed)
    X_scaled = scaler.transform(X_transformed)

    kf_inner = KFold(n_splits=10, shuffle=True, random_state=0)
    kf_outer = KFold(n_splits=10, shuffle=True, random_state=1)

    print "{} Begin Model Fitting and cross-validation for individual feature: {}".format(timestamp(), feature_name)
    # Set the parameters by cross-validation:
    tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4], 'C': [0.001, 0.01, 0.1, 1, 10, 100, 1000]}]
    clf = GridSearchCV(SVR(), tuned_parameters, cv=kf_inner, scoring=scoreType, n_jobs=10)
    clf.fit(X_scaled, y)
    print "{} Finished GridSearch and fit for individual feature: {}".format(timestamp(), feature_name)
    print "{} Best cross-validated mean score of the best_estimator during grid search: {:.3f}".format(timestamp(), clf.best_score_)
    print "{} Mean RMSE: {:.1f} C".format(timestamp(), np.sqrt(-1*clf.best_score_))
    
    scores_best_est = cross_val_score(clf, X_scaled, y, cv=kf_outer, scoring=scoreType, n_jobs=10)
    print "{} Finished fit and scoring for individual feature:    {}".format(timestamp(), feature_name)
    print "{} Mean cross-validated score of the best_estimator during cross-val scoring: {:.3f}".format(timestamp(), np.mean(scores_best_est))
    print "{} Mean RMSE: {:.1f} C".format(timestamp(), np.sqrt(-1*np.mean(scores_best_est))) 
    
    for score in scores_best_est:
        outFileHandle.write(feature_name+','+str(score)+'\n')
        outFileHandle.flush()                          
    data_Dict[feature_name] = scores_best_est

    # change to cross_val_predict
    y_pred = cross_val_predict(clf, X_scaled, y, cv=kf_outer, n_jobs=10)
    y_true = list(y)
    
    return data_Dict, y_true, y_pred

# Now run the algorithm and get the scores/predictions:  

outFileHandle = open(outputFile, 'w')
outFileHandle.write("FEATURE,SCORE"+"\n")

# Get a score for each feature individually (except dipeptide frequencies):
query = "select * from GOLD_FEATURES where OGT is not null;"
c.execute(query)
results = c.fetchall()
dataArray = np.asarray(results)
headerNames = [description[0] for description in c.description]
non_features = set(['GOLD_ID', 'OGT', 'NORMALIZED_DISULFIDE_BONDS', 'F_DISULFIDE_BONDS', 'HTPP', 'NORMALIZED_SALT_BRIDGES', 'OGT_RANGE', 'OGT_BIN', 'rRNA_16S_GC'])
non_features_indices = set([])
for i, feature in enumerate(headerNames):
    if feature in non_features:
        non_features_indices.add(i)

dataDict = {}
for index in range(len(headerNames)): 
    if index not in non_features_indices:
        X_feature = dataArray[:, index:index+1].astype(float) 
        y_feature = dataArray[:, 1].astype(float)
        featureLabel = headerNames[index]
        dataDict, y_true_feature, y_pred_feature = FeaturePrediction(X_feature, y_feature , featureLabel, dataDict)

# Then for each feature combined, but not including dipeptide:
dataDict = {}
idx_in_columns = [i for i in range(len(headerNames)) if i not in non_features_indices]
X_combined = dataArray[:, idx_in_columns].astype(float)
y_combined = dataArray[:, 1].astype(float)
dataDict, y_true_combined, y_pred_combined = FeaturePrediction(X_combined, y_combined, 'COMBINED_WO_DIPEPTIDE', dataDict)

for score in dataDict['COMBINED_WO_DIPEPTIDE']:
    outFileHandle.write('COMBINED_WO_DIPEPTIDE,'+str(score)+'\n')


# Then for dipeptide frequencies as a whole:
dataDict = {}
query = "select * from DIPEPTIDE_FREQUENCIES where OGT is not null;"
c.execute(query)
results = c.fetchall()
dataArray = np.asarray(results)
X_dipep = dataArray[:, 2:].astype(float) # Exclude column 0 (GOLD_ID) and column 1 (OGT)
y_dipep = dataArray[:, 1].astype(float) # Get column of OGTs (column 1)
dataDict, y_true_dipep, y_pred_dipep = FeaturePrediction(X_dipep, y_dipep, 'DIPEPTIDE_FREQUENCIES', dataDict)

for feature in dataDict:
    for score in dataDict[feature]:
        outFileHandle.write(feature+','+str(score)+'\n')

# Write predictions to file:
predHandle_dipep = open(predFile_dipep, 'w')
predHandle_dipep.write('ACTUAL'+','+'PREDICTED'+'\n')
for i in range(len(y_true_dipep)):
    predHandle_dipep.write(str(y_true_dipep[i])+','+str(y_pred_dipep[i])+'\n')

predHandle_dipep.close()

predHandle_combined = open(predFile_combined_wo_dipep, 'w')
predHandle_combined.write('ACTUAL'+','+'PREDICTED'+'\n')
for i in range(len(y_true_combined)):
    predHandle_combined.write(str(y_true_combined[i])+','+str(y_pred_combined[i])+'\n')

predHandle_combined.close()
outFileHandle.close()
conn.close() 

