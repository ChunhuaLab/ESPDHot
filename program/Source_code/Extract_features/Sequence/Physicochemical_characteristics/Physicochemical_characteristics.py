# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 15:33:51 2022

@author:Lianci Tao
"""
import pandas as pd
import numpy as np
import sys

# input PDB file


file = sys.argv[1].split('/')[-1]
pdbID = file.split('.')[0]
chain = pdbID.split('_')[1]
pdb = open(sys.argv[1], 'r')
content = pdb.read()
pdb.close()
line = content.strip().split("\n")
# table lookup
ATOM = []
F_i = []
for i in range(0,len(line)):
    if line[i][0:4] == 'ATOM' and line[i][21] == chain and ((line[i][13] == 'C' and line[i][14] == 'A')):
        ATOM_split = line[i].split()
        ATOM.append(ATOM_split)
        pdbname = line[i][17:20]
        if pdbname =='GLY':
            resname = 'G'
            F = [5.97,57,-1.62]
        elif pdbname =='ALA':
            resname = 'A'
            F = [6,71.1,	-0.22]
        elif pdbname =='VAL':
            resname = 'V'
            F = [5.96,99.1,4.45]
        elif pdbname == 'LEU':
            resname = 'L'
            F = [5.98,113.2,5.01]
        elif pdbname == 'ILE':
            resname = 'I';
            F = [6.02,113.2,5.58]
        elif pdbname == 'PRO':
            resname = 'P'
            F = [6.3,97.1,-3.03]
        elif pdbname == 'PHE':
            resname = 'F'
            F = [5.48,147.2,5.27]
        elif pdbname == 'TRP':
            resname = 'W';
            F = [5.89,186.2,5.2]
        elif pdbname == 'TYR':
            resname = 'Y'
            F = [5.66,163.2,2.15]
        elif pdbname == 'SER':
            resname = 'S'
            F = [5.68,87.1,-2.84]
        elif pdbname == 'THR':
            resname = 'T'
            F = [5.66,101.1,-1.2]
        elif pdbname == 'CYS':
            resname = 'C'
            F = [5.05,103.1,4.66]
        elif pdbname == 'MET':
            resname = 'M'
            F = [5.74,131.2,3.51]
        elif pdbname == 'ASN':
            resname = 'N'
            F = [5.41,114.1,-2.65]
        elif pdbname == 'GLN':
            resname = 'Q'
            F = [5.65,128.1,-2.76]
        elif pdbname == 'ASP':
            resname = 'D'
            F = [2.77,115.1,-4.12]
        elif pdbname == 'GLU':
            resname = 'E'
            F = [3.22,129.1,-3.64]
        elif pdbname == 'HIS':
            resname = 'H'
            F = [7.59,137.1,1.28]
        elif pdbname == 'LYS':
            resname = 'K'
            F = [9.74,128.2,-4.18]
        elif pdbname == 'ARG':
            resname = 'R'
            F = [10.76,156.2,-0.93]
        F_i.append(F)
F_i = np.array(F_i)
pos = [[x[5]] for x in ATOM]
pos = np.array(pos)
data = pd.DataFrame(np.hstack((pos,F_i)))
data.columns = ['pos','Isoep','Mass','Enc']
data.to_csv(pdbID+'_PC.csv',index=False)
