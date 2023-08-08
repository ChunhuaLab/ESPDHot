# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 15:33:51 2022

@author:Lianci Tao
"""
import sys,re
import pandas as pd
import numpy as np


file = sys.argv[1].split('/')[-1]
pdbID = file.split('.')[0]
pdbnum = pdbID.split('_')[0]
chain = pdbID.split('_')[1]
data = []
with open(sys.argv[1], 'r') as f:
    content = f.readlines()
    f.close()

data = pd.DataFrame(np.array([re.split(r'\s+', ct_i.strip()) for ct_i in content]))

out_name = r'./' +pdbID+'.csv'
data.to_csv(out_name,header=False,index=False)







