# -*- coding: utf-8 -*-
"""
Created on Nov 5, 2022

@author: Lianci Tao
"""

# import resources
import sys
import os
import pandas as pd
import numpy as np
import joblib


# data normalization
def testdataPrepare(X, mean, std):
    nor_X = (X - mean) / std
    return nor_X


# model performance
def performance(labelArr, predictArr):
    TP = 0.
    TN = 0.
    FP = 0.
    FN = 0.
    for i in range(len(labelArr)):
        if labelArr[i] == 1 and predictArr[i] == 1:
            TP += 1.
        if labelArr[i] == 1 and predictArr[i] == 0:
            FN += 1.
        if labelArr[i] == 0 and predictArr[i] == 1:
            FP += 1.
        if labelArr[i] == 0 and predictArr[i] == 0:
            TN += 1.
    ACC = (TP + TN) / (TP + FN + FP + TN)
    SEN = TP / (TP + FN + 0.0000001)
    SPE = TN / (FP + TN + 0.0000001)
    PRE = TP / (TP + FP + 0.0000001)
    F1 = 2 * (PRE * SEN) / (PRE + SEN + 0.0000001)
    fz = float(TP * TN - FP * FN + 0.0000001)
    fm = float(np.math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) + 0.0000001))
    MCC = fz / fm
    return ACC, SEN, SPE, PRE, F1, MCC


if __name__ == "__main__":

    # file
    file = sys.argv[1].split('/')[-1]
    pdbID = file.split('_f')[0]
    outdir = "./Result/"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # input mean and std for data normalization
    mean = joblib.load('./Model_file/mean.param')
    std = joblib.load('./Model_file/std.param')

    # input RSSHOT method
    eclf = joblib.load('./Model_file/model.m')

    # input test data
    testing_dataset = pd.read_csv(sys.argv[1], encoding='utf-8', sep=',', quoting=1)
    data = testing_dataset.drop(columns=['pos', 'res'])
    tmp = testing_dataset[['pos', 'res']]

    # data normalization
    prepared_testset = testdataPrepare(data, mean, std)
    X = np.array(prepared_testset)

    # prediction
    y_pred = eclf.predict(X)
    proba = eclf.predict_proba(X)
    pre_prob = []
    for i in range(len(proba)):
        pre_prob.append(proba[i][1])
    pre_prob = [round(i, 3) for i in pre_prob]

    tmp = np.array(tmp)
    pre_prob = np.array(pre_prob).reshape(-1, 1)
    y_pred = y_pred.reshape(-1, 1)
    result = pd.DataFrame(np.hstack((tmp, y_pred, pre_prob)))
    result.columns = ['Position', 'Residue', 'Predict_label', 'Predict_probability']
    result.to_csv(outdir + pdbID + '_hotspot_result.csv', index=False)

    print('\nThe prediction is done.\n')
    print('The result is in the file "' + outdir + pdbID + '_hotspot_result.csv".')

