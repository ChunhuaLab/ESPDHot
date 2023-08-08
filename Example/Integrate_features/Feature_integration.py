# -*- coding: utf-8 -*-
"""
Created on Monday nov  14 15:33:51 2022

@author: Lianci Tao
"""

import pandas as pd
import sys
import os
import numpy as np
import re
# input PDB file
outdir = "./Output/"
if not os.path.exists(outdir):
    os.makedirs(outdir)
file = sys.argv[1].split('/')[-1]
pdbID = file.split('.')[0]
pdbnum = pdbID.split('_')[0]
chain = pdbID.split('_')[1]
pdb = open(sys.argv[1], 'r')
content = pdb.read()
pdb.close

line = content.strip().split("\n")
# confirm number of 1st residue in FASTA file

for i in range(0, len(line)):
    str1 = re.sub(" +", ' ', str(line[i])).split(' ')
    if str1[0] == 'ATOM' and str1[2] == 'CA' and str1[4] == chain:
        Seq1stNum = int(str1[5])
        break

# Physicochemical characteristics
PC = pd.read_csv("./Input/" + pdbID + "_PC.csv",
                 encoding='utf-8', sep=',', quoting=1)

# PSSM
PSSM = pd.read_table("./Input/" + pdbID + ".asn_matrix.txt",
                     encoding='utf-8', sep='\t+', quoting=1, engine='python',)

PSSM['P'] += (Seq1stNum-1)
PSSM = PSSM[['P', 'Q', 'M', 'N', 'R', 'G']]
PSSM.columns = ['pos', 'PQ', 'PM', 'PN', 'PR', 'PG']
position = PSSM['pos']

##coevolution
coevo = pd.read_csv('./Input/'+pdbID+"_coevolution.csv", encoding = 'utf-8', sep = ',', quoting = 1)
coevo['pos'] = position
# CX/DPX
tbl_bound = pd.read_table("./Input/" + pdbnum + "_bound.tbl",
                          encoding='utf-8', sep=' +', quoting=1,
                          engine='python', skiprows=9, header=None)
tbl_bound.columns = ['chain', 'pos', 'AA',
                     'BmDPXa', 'BsdDPXa', 'BmDPXs', 'BsdDPXs', 'BmaxDPX', 'BminDPX',
                     'BmCXa', 'BsdCXa', 'BmCXs', 'BsdCXs', 'BmaxCX', 'BminCX']
tbl_bound = tbl_bound[['chain', 'pos', 'AA', 'BsdDPXa', 'BmDPXs', 'BmCXa', 'BsdCXa']]
tbl_bound = tbl_bound.loc[tbl_bound['chain'] == chain]

tbl_unbound = pd.read_table("./Input/" + pdbnum + "_unbound.tbl",
                            encoding='utf-8', sep=' +', quoting=1,
                            engine='python', skiprows=9, header=None)
tbl_unbound.columns = ['chain', 'pos', 'AA',
                       'UmDPXa', 'UsdDPXa', 'UmDPXs', 'UsdDPXs', 'UmaxDPX', 'UminDPX',
                       'UmCXa', 'UsdCXa', 'UmCXs', 'UsdCXs', 'UmaxCX', 'UminCX']
tbl_unbound = tbl_unbound.loc[tbl_unbound['chain'] == chain]
tbl_unbound = tbl_unbound[['pos', 'UsdDPXa', 'UmDPXs', 'UsdCXa']]

tbl = pd.merge(tbl_bound, tbl_unbound, on=['pos'])

tbl['DsdDPXa'] = tbl['BsdDPXa'] - tbl['UsdDPXa']
tbl['DmDPXs'] = tbl['BmDPXs'] - tbl['UmDPXs']
tbl['DsdCXa'] = tbl['BsdCXa'] - tbl['UsdCXa']
tbl = tbl[['pos', 'DsdDPXa', 'DmDPXs', 'BmCXa', 'DsdCXa']]

##secondary structure
ss = pd.read_excel("./Input/" + pdbID + "_structure.xls")

ss.columns=['seq','res','structure','ASA','Phi','Psi','Theta','Tau', 'HSE_A_UP', 'HSE_A_down','probC', 'probH', 'probE']
ss = ss[['Tau', 'Theta','probC', 'probH', 'probE']]
ss['pos'] = position

# ASA
ASA_bound = open("./Input/" + pdbID + "_bound.rsa", 'r')
content = ASA_bound.read()
ASA_bound.close()
line1 = content.strip().split("\n")
ASA_bound = []
for i in range(0, len(line1)):
    if line1[i][0:3] == 'RES' and line1[i][8] == chain:
        ASA_bound_split = line1[i].split()
        ASA_bound.append(ASA_bound_split)
ASA_bound = [[x[3], x[10]] for x in ASA_bound]
ASA_bound = pd.DataFrame(ASA_bound)
ASA_bound.columns = ['pos', 'BaSASAnp']

ASA_unbound = open("./Input/" + pdbID + '_unbound.rsa', 'r')
content = ASA_unbound.read()
ASA_unbound.close()
line2 = content.strip().split("\n")
ASA_unbound = []
for i in range(0, len(line2)):
    if line2[i][0:3] == 'RES' and line2[i][8] == chain:
        ASA_unbound_split = line2[i].split()
        ASA_unbound.append(ASA_unbound_split)
ASA_unbound = [[x[3],x[10]] for x in ASA_unbound]
ASA_unbound = pd.DataFrame(ASA_unbound)
ASA_unbound.columns = ['pos', 'UaSASAnp']

ASA = pd.DataFrame(ASA_unbound['UaSASAnp'].apply(float)- ASA_bound['BaSASAnp'].apply(float))
ASA.columns=['DaSASAnp']
ASA = pd.concat([ASA_bound['pos'], ASA], axis=1)
ASA['pos'] = ASA['pos'].apply(int)

# Solvent exposure
HSE = pd.read_table("./Input/" + pdbID + "_HSE.txt",
                    encoding='utf-8', quoting=1, sep=' +', engine='python')
HSE['Position'] = range(Seq1stNum, len(HSE) + Seq1stNum)
HSE.columns = ['pos', 'res', 'HSEup', 'HSEdown', 'CN']

# IP
Interface_P = pd.read_csv("./Input/" + pdbID + "_IP.csv",
                 encoding='utf-8', sep=',', quoting=1)
Interface_P =pd.DataFrame(Interface_P )
# ANN
ANN_H = pd.read_csv("./Input/" + pdbID + "_Hydrophobicity.txt",
                    encoding='utf-8', sep=' ', quoting=1)
ANN_H = ANN_H[['Resid', 'Bw', 'Cw']][ANN_H['chain'] == chain]
ANN_H.columns = ['pos', 'Bh', 'Ch']
ANN_H = ANN_H.reset_index(drop=True)
ANN_S = pd.read_csv("./Input/" + pdbID + "_SAS.txt",
                    encoding='utf-8', sep=' ', quoting=1)
ANN_S = ANN_S[['Resid', 'Kw']][ANN_S['chain'] == chain]
ANN_S.columns = ['pos', 'Ks']
ANN_S = ANN_S.reset_index(drop=True)

ANN_P = pd.read_csv("./Input/" + pdbID + "_Polarity.txt",
                    encoding='utf-8', sep=' ', quoting=1)

ANN_P = ANN_P[['Resid', 'Cw']][ANN_P['chain'] == chain]
ANN_P.columns = ['pos', 'Cp']
ANN_P = ANN_P.reset_index(drop=True)
ANN_M = pd.read_csv("./Input/" + pdbID + "_Mass.txt",
                    encoding='utf-8', sep=' ', quoting=1)
ANN_M = ANN_M[['Resid', 'Cw']][ANN_M['chain'] == chain]
ANN_M.columns = ['pos', 'Cm']
ANN_M = ANN_M.reset_index(drop=True)
ANN1 = pd.merge(ANN_M, ANN_H,on = ['pos'])
ANN2 = pd.merge(ANN1, ANN_S,on = ['pos'])
ANN = pd.merge(ANN2, ANN_P,on = ['pos'])
##Dynamics
Dynamics = pd.read_csv("./Input/" + pdbID + "_dynamics.csv",
                       encoding='utf-8', sep=',', quoting=1)
Dynamics = Dynamics[['pos', 'F4', 'K3', 'F5', 'K10']]

# Integration
M1 = pd.merge(PC, PSSM, on=['pos'])

M2 = pd.merge(M1, coevo, on=['pos'])

M3 = pd.merge(M2, tbl, on=['pos'])

M4 = pd.merge(M3, ss, on=['pos'])

M5 = pd.merge(M4, ASA, on=['pos'])

M6 = pd.merge(M5, HSE, on=['pos'])
M7 = pd.merge(M6, ANN, on=['pos'])
M8 = pd.merge(M7, Dynamics, on=['pos'])

M9 = pd.merge(M8, Interface_P, on=['pos','res'])

res = M9[['pos', 'res', 'Isoep', 'Mass', 'PM', 'PQ', 'PR', 'Theta', 'probC', 'probH',
          'probE', 'DsdDPXa', 'DmDPXs', 'HSEup', 'CN', 'DaSASAnp', 'MI', 'DI', 'Cm',
          'Bh', 'Ch', 'Ks', 'F4', 'K3', 'Enc', 'IP', 'Tau', 'BmCXa', 'DsdCXa', 'HSEdown',
          'Cp', 'F5', 'PN', 'K10', 'PG']]

res.to_csv(outdir + pdbID + "_feature.csv", index=False)

