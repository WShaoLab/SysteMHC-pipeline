#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

print ('usage: please input argv "argv1:peptidelist", "argv2:ionslist","argv3:MHCpan.xls","argv4:outname"')

pep = pd.read_csv(sys.argv[1],header=None, sep = '\t',engine='python')
pep.columns = ['Peptide']
ions = pd.read_csv(sys.argv[2],header=None, sep = '\t',engine='python')
ions.columns = ['ions']

damhc = pd.read_csv(sys.argv[3], header=None, sep = '\t',engine='python')
damhc.columns = damhc.iloc[1,:]
if 'Pos' in damhc.index:
    damhc.index= np.arange(len(damhc.index))

da = damhc.drop([0,1],axis=0)
da.NB = pd.to_numeric(da.NB)
bind = da[da.NB>0].loc[:,['Peptide','NB']]

data = pd.concat([pep,ions],axis=1)
data1 = pd.merge(bind,data,on='Peptide')
bindions = data1.loc[:,'ions']
bindpeps = bind.loc[:,'Peptide']
outname1 = sys.argv[4] + '_bindions.csv'
outname2 = sys.argv[4] + '_bindpeps.csv'
bindions.to_csv(outname1,header=None,index=False)
bindpeps.to_csv(outname2,header=None,index=False)

