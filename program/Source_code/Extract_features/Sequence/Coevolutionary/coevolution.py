# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 15:33:51 2022

@author:Lianci Tao
"""
import pandas as pd
import numpy as np
import sys

# input coecolutionary file
file = sys.argv[1].split('/')[-1]
pdbID = file.split('.')[0]
coe = open(sys.argv[1], 'r')
content = coe.read()
coe.close()
line = content.strip().split("\n")
# table lookup
nlen = int(line[len(line)-1].split()[0])+1
MI =np.zeros(nlen,dtype=float)
DI =np.zeros(nlen,dtype=float)
for num_i in np.arange(0,nlen):
    for i in range(0,len(line)):
        if num_i == int(line[i].split()[0])-1 or num_i == int(line[i].split()[1])-1:

            MI[num_i] += float(line[i].split()[2])
            DI[num_i] += float(line[i].split()[3])
sum_MI = sum(MI)
sum_DI = sum(DI)
MI = np.array(MI)/(sum_MI/float(nlen))
DI = np.array(DI)/(sum_DI/float(nlen))
coevolutonary = pd.concat([pd.DataFrame(MI),pd.DataFrame(DI)],axis=1)
coevolutonary.columns=['MI','DI']
print(pdbID)
coevolutonary.to_csv(pdbID+'_coevolution.csv',index=False)


















