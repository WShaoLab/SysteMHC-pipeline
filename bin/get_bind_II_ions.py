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

damhc0 = damhc.iloc[0:2,:]
damhc1 = damhc.iloc[2:,:]
indx0 = damhc0.index
indx1 = damhc1.index
a0 = list(indx0[1])
a1 = []
for i in range(damhc1.shape[0]):
    ax = list(indx1[i])
    a1.append(ax)

da1 = pd.DataFrame(a1, columns=a0)
da2 = damhc1.copy(deep=True)
da2.index = np.arange(damhc1.shape[0])
da2.columns = list(damhc0.iloc[1,:])
da = pd.concat([da1,da2], axis=1)

da.NB = pd.to_numeric(da.NB)
bindx = da[da.NB>0]
bind = bindx.loc[:,['Peptide','NB']]

data = pd.concat([pep,ions],axis=1)
data1 = pd.merge(bind,data,on='Peptide')
bindions = data1.loc[:,'ions']
bindpeps = bind.loc[:,'Peptide']
outname1 = sys.argv[4] + '_bindions.csv'
outname2 = sys.argv[4] + '_bindpeps.csv'
bindions.to_csv(outname1,header=None,index=False)
bindpeps.to_csv(outname2,header=None,index=False)

