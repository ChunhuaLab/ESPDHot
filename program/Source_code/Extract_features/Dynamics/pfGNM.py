# -*- coding: utf-8 -*-
"""
Created on 11.14 2022
Attention:this is pfGNM

@author: Tong Zhou

"""

import numpy as np
import pandas as pd
import sys


def readpdb(line):
    ATOM = []
    for i in range(0, len(line)):
        if line[i][0:4] == 'ATOM' and line[i][21] == chain and ((line[i][13] == 'C' and line[i][14] == 'A')):
            ATOM_split = line[i].split()
            ATOM.append(ATOM_split)
    coordinate = [[eval(x[6]), eval(x[7]), eval(x[8])] for x in ATOM]
    pos = [[x[5]] for x in ATOM]
    return coordinate, pos


def distance(x1, x2):
    dis = np.sqrt(np.square(x1[0] - x2[0]) + np.square(x1[1] - x2[1]) + np.square(x1[2] - x2[2]))
    return dis


def mainfuc():
    coordinate, pos = readpdb(line)
    n = len(coordinate)
    netmat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            else:
                dis = distance(coordinate[i], coordinate[j])
                x = 1 / np.square(dis)
                netmat[i, j] = -x
                netmat[j, i] = -x

    for i in range(n):
        netmat[i, i] = -1 * sum(netmat[i])

    fluslow1 = []
    D, V = np.linalg.eig(netmat)
    D = np.array(D)
    V = np.array(V)
    sorted_indices = np.argsort(D)
    D=D[sorted_indices]
    V =V[:, sorted_indices]

    for i in range(n):
        tmp = V[i, 1] * V[i, 1] / D[1]
        fluslow1.append(tmp)

    f1 = np.array(np.around(fluslow1, 4)).reshape(-1, 1)
    f1max = max(f1)
    f1 = f1/f1max

    fluslow2 = []
    for i in range(n):
        tmp = 0
        for j in range(1, 3):
            tmp = tmp + V[i, j] * V[i, j] / D[j]
        fluslow2.append(tmp)
    f2 = np.array(np.around(fluslow2, 4)).reshape(-1, 1)
    f2max = max(f2)
    f2 = f2 / f2max

    fluslow3 = []
    for i in range(n):
        tmp = 0
        for j in range(1, 4):
            tmp = tmp + V[i, j] * V[i, j] / D[j]
        fluslow3.append(tmp)
    f3 = np.array(np.around(fluslow3, 4)).reshape(-1, 1)
    f3max = max(f3)
    f3 = f3 / f3max

    fluslow4 = []
    for i in range(n):
        tmp = 0
        for j in range(1, 5):
            tmp = tmp + V[i, j] * V[i, j] / D[j]
        fluslow4.append(tmp)
    f4 = np.array(np.around(fluslow4, 4)).reshape(-1, 1)
    f4max = max(f4)
    f4 = f4 / f4max

    fluslow5 = []
    for i in range(n):
        tmp = 0
        for j in range(1, 6):
            tmp = tmp + V[i, j] * V[i, j] / D[j]
        fluslow5.append(tmp)
    f5 = np.array(np.around(fluslow5, 4)).reshape(-1, 1)
    f5max = max(f5)
    f5 = f5 / f5max

    fluslow6 = []
    for i in range(n):
        tmp = 0
        for j in range(1, 7):
            tmp = tmp + V[i, j] * V[i, j] / D[j]
        fluslow6.append(tmp)
    f6 = np.array(np.around(fluslow6, 4)).reshape(-1, 1)
    f6max = max(f6)
    f6 = f6 / f6max

    flufast1 = []
    for i in range(n):
        tmp = 0
        for j in range(n-1, n):
            tmp = tmp + V[i, j] * V[i, j] / D[j]
        flufast1.append(tmp)
    K1 = np.array(np.around(flufast1, 30)).reshape(-1, 1)
    K1max = max(K1)
    K1 = K1 / K1max

    flufast3 = []
    for i in range(n):
        tmp = 0
        for j in range(n-3, n):
            tmp = tmp + V[i, j] * V[i, j] / D[j]
        flufast3.append(tmp)
    K3 = np.array(np.around(flufast3, 30)).reshape(-1, 1)
    K3max = max(K3)
    K3 = K3 / K3max

    flufast10 = []
    for i in range(n):
        tmp = 0
        for j in range(n - 10, n):
            tmp = tmp + V[i, j] * V[i, j] / D[j]
        flufast10.append(tmp)
    K10 = np.array(np.around(flufast10, 30)).reshape(-1, 1)
    K10max = max(K10)
    K10 = K10 / K10max



    dynamics = pd.DataFrame(np.hstack((pos, f1, f2, f3, f4, f5, f6,K1,K3,K10)))
    dynamics.columns = ['pos', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6','K1','K3','K10']
    return dynamics


if __name__ == "__main__":

    file = sys.argv[1].split('/')[-1]
    pdbID = file.split('.')[0]
    chain = pdbID.split('_')[1]
    pdb = open(sys.argv[1], 'r')
    content = pdb.read()
    pdb.close()
    line = content.strip().split("\n")
    dynamics = mainfuc()
    dynamics.to_csv(pdbID + '_dynamics.csv', index=False)



