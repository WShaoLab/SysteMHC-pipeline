#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

#print ('usage: please input argv "argv1: merged.pep.xml converted merged.csv"')

def cal_fdr(inputda):
    data = inputda.copy()
    data['fdr']=0
    data.index=range(data.shape[0])
    ne=0
    for i in range(data.shape[0]):
        countn = i+1
        if 'DECOY_' in data.loc[i,'Protein'] :
            ne += 1
        fdr = ne/countn
        data.loc[i,'fdr'] = np.round(fdr,4)
    return data

def getprob(PSM,peptide):
    
    for i in range(PSM.shape[0]-1,0,-1):
        if (PSM.loc[i,'fdr'] <= 0.01):
            pb1 = PSM.loc[i,'Probability']
            break

    for i in range(peptide.shape[0]-1,0,-1):
        if (peptide.loc[i,'fdr'] <= 0.01):
            pb2 = peptide.loc[i,'Probability']
            break
    pb = max(pb1,pb2)
    return pb1,pb2,pb

data = pd.read_csv(sys.argv[1],sep=',')
data = data.sort_values(by=['Probability'],ascending=False)
PSM = cal_fdr(data)
data1 = data.drop_duplicates(subset=['Peptide_Sequence'],keep='first')
peptide = cal_fdr(data1)

pb1,pb2,pb = getprob(PSM,peptide)
data2 = PSM.loc[PSM.Probability >= pb]
data3 = peptide.loc[peptide.Probability >= pb]
data4 = data3.loc[(data3['Peptide_Sequence'].str.len()>=8) & (data3['Peptide_Sequence'].str.len()<=14)]
peplist = data4.loc[:,'Peptide_Sequence']
ionslist = data4.loc[:,'Ions']

data2.to_csv('psm_result.csv',index=False)
data3.to_csv('peptide_result.csv',index=False)
data4.to_csv('peptide_814_result.csv',index=False)
peplist.to_csv('peplist.csv',header=None,index=False) 
ionslist.to_csv('ionslist.csv',header=None,index=False) 

print(pb)
