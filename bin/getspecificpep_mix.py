#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt

#print ('usage: please input argv "argv1:peptidelist", "argv2:ionslist","argv3:MHCpan.xls","argv4:outname"')

def readc2(file1,alle):    
    df = []
    df2 = []
    
    for i in range(len(file1)):
        dat  = []
        inp1 = os.path.join(dirname,file1[i])

        da1 = pd.read_csv(inp1,sep='\t')
        da2 = da1.copy()
        da2['NB'] = np.nan
        da2.loc[da2['%Rank_best']<=5,['NB']]=1
        da3 = da2[da2['NB']==1]
        col1 = ['Peptide','%Rank_best','NB']
        da4 = da3[col1]
        da4.columns = ['Peptide',alle[i],'NB']
        df.append(da4)
            
    dax = df[0]
    
    if len(alle)>=2 :
        for i in range(1,len(alle)):
            dax = pd.merge(dax,df[i],on='Peptide')
    return dax
    
def calculate(dax,alle,dap):
    
    #dalx['annotation_score'] = round(x.max(axis=1)/x.min(axis=1)).astype(int)

    col = ['Peptide']
    for i in dax.columns:
        if 'NB' in i:
            col.append(i)

    dax1 =dax[(dax==1).any(axis=1)]

    daf=dax1.drop(columns=col)

    daf2 = daf.astype(float)
    topa = daf2.idxmin(axis = 1)

    dax2 = dax1.copy(deep=True)
    dax2['top_allele'] = topa

    dax3 = dax2.loc[:,['Peptide','top_allele']]
    dax3.columns = dax3.columns.str.replace('Peptide','Peptide_Sequence')
    
    dan = pd.merge(dap,dax3,how='outer',on=['Peptide_Sequence'])
    dan['top_allele']= dan.loc[:,'top_allele'].fillna('unclassified')
    dan['annotation_score']=np.nan

    return dan


path2 = sys.argv[1]
path1 = os.getcwd()
dirname = os.path.join(path1,path2)

#filepath2 = 'Results/MHCSpecificLib/netmhcout/'
#dirname = res + filepath2

#filepath1 = 'Results/Peplevel_FDR/'
#infilepath = res + filepath1
#infile = infilepath + 'peptide_result.csv'
infile = sys.argv[2]

if os.path.exists(dirname):
    
    anno_file = 'alle_annotation.csv'
    
    if not os.path.exists(anno_file):
       
        files = os.listdir(dirname)
        file=[]
        alle=[]
        file2 = []

        for i in files:
            if 'csv' in i:
                file.append(i)
        file1 = sorted(file)

        for i in file1:

            dirStr, ext = os.path.splitext(i)
            alle.append(dirStr)
             
        dapx = pd.read_csv(infile)
        dap = dapx.loc[(dapx['Peptide_Sequence'].str.len()>=9) & (dapx['Peptide_Sequence'].str.len()<=25)]

        dax = readc2(file1,alle)
        danx = calculate(dax,alle,dap)
        dan = danx.copy()
        dan['length'] = dan['Peptide_Sequence'].str.len()

        out3 = 'classII_peptide_complete.csv'
        dan.to_csv(out3,index=False)

