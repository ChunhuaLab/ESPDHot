# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 15:33:51 2022

@author:Lianci Tao
"""
import pandas as pd
import numpy as np
import sys
tmp_res_i = []
pdb_IP = []
res_i=[]
file = sys.argv[1].split('/')[-1]
pdbID = file.split('.')[0]

dssp = open(sys.argv[1], 'r', encoding='UTF-8')
content = dssp.read()
dssp.close()
Dssp = content.strip().split("\n")##将dssm文本按行读取

for jk in np.arange(0,len(Dssp)):
    tmp = Dssp[jk].split()#每行按空格进行分割
    if tmp[0]=='#':
        num = jk+1 #找到dssp有用数据起始行

for jc in np.arange(num,len(Dssp)):
    tmp = Dssp[jc].split()
    res_i.append(tmp[1])
    tmp_res = tmp[3]

    tmp_struct = tmp[4]

    if tmp_struct == 'B'or tmp_struct=='G'or tmp_struct=='T':#把八种结构分为X、Y、Z
        stru = 'X'
        res_P = open(sys.argv[2], 'r',encoding='UTF-8')
        content_p = res_P.read()
        res_P.close()
        RES_P = content_p.strip().split("\n")  ##将dssp文本按行读取
        for i in np.arange(1, len(RES_P)):
            find_P = RES_P[i].split()# 每行按空格进行分割

            if find_P[0] == tmp_res:
                if find_P[1]== stru:

                    IP = find_P[2]
                    tmp_res_i.append(tmp_res)
                    pdb_IP.append(IP)


    elif tmp_struct == 'E'or'I':
        stru = 'Z'
        res_P = open(sys.argv[2], 'r',encoding='UTF-8')
        content_p = res_P.read()
        res_P.close()
        RES_P = content_p.strip().split("\n")  ##将dssp文本按行读取
        for i in np.arange(1, len(RES_P)):
            find_P = RES_P[i].split()# 每行按空格进行分割
            if find_P[0] == tmp_res:
                if find_P[1]== stru:
                    IP = find_P[2]
                    tmp_res_i.append(tmp_res)
                    pdb_IP.append(IP)
    else:
        stru = 'Y'
        res_P = open(sys.argv[2], 'r',encoding='UTF-8')
        content_p = res_P.read()
        res_P.close()
        RES_P = content_p.strip().split("\n")  ##将dssp文本按行读取
        for i in np.arange(1, len(RES_P)):
            find_P = RES_P[i].split()# 每行按空格进行分割
            if find_P[0] == tmp_res:
                if find_P[1]== stru:
                    IP = find_P[2]
                    tmp_res_i.append(tmp_res)
                    pdb_IP.append(IP)

tmp_res_i = np.array(tmp_res_i)
pdb_IP = np.array(pdb_IP)
res_i= np.array(res_i)
data = pd.DataFrame(np.transpose(np.vstack((res_i,tmp_res_i,pdb_IP))))
data.columns = ['pos','res','IP']
data.to_csv(pdbID+'_IP.csv',index=False)

