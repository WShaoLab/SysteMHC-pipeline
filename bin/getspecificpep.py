#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt

#print ('usage: please input argv "argv1:peptidelist", "argv2:ionslist","argv3:MHCpan.xls","argv4:outname"')

def readc1(file1,alle,file2):
    df =[]
    df2 = []

    for i in range(len(file1)):

        inp1 = os.path.join(dirname,file1[i])
        damhc = pd.read_csv(inp1, header=None, sep = '\t',engine='python')
        damhc.columns = damhc.iloc[1,:]
        if 'Pos' in damhc.index:
            damhc.index= np.arange(len(damhc.index))
        da = damhc.drop([0,1],axis=0)
        da.NB = pd.to_numeric(da.NB)
        da1 = da.loc[:,['Peptide','EL_Rank','NB']]
        da1.columns = ['Peptide',alle[i],'NB']
        df.append(da1)
    
        linea=[]
        linea1=[]

        inp2 = os.path.join(dirname,file2[i])
        with open(inp2,'r') as f:
            lines = f.read().splitlines()
            for line in lines:
                if (re.match(r'^\s\s\s\d',line)):
                    linea.append(line)
                    nm = line.split(' ')
                    ax = list(filter(None,nm))
                    linea1.append(ax)
                
            dal = pd.DataFrame(linea1)  
            dal1 = dal.iloc[:,[2,15]]
            dal1.columns = ['Peptide',alle[i]]
            df2.append(dal1)
        
    dax = df[0]
    dalx = df2[0]
    if len(alle)>=2 :
        for i in range(1,len(alle)):
            dax = pd.merge(dax,df[i],on='Peptide')
            dalx = pd.merge(dalx,df2[i],on='Peptide')
            
    return dax,dalx 

def readc2(file1,alle,file2):    
    df = []
    df2 = []
    
    for i in range(len(file1)):
        dat  = []
        inp1 = os.path.join(dirname,file1[i])

        with open(inp1, 'r',encoding='utf-8-sig') as f_input:
            for line in f_input:
                dat.append(list(line.strip().split('\t')))
        damhc = pd.DataFrame(dat)
        damhc0 = damhc.iloc[1,:]
        damhc1 = damhc.iloc[2:,:]
        damhc1.columns = damhc0
        damhc1.index = np.arange(damhc1.shape[0])
        da = damhc1.copy(deep=True)
        da.NB = pd.to_numeric(da.NB)
        da1 = da.loc[:,['Peptide','Rank','NB']]
        da1.columns = ['Peptide',alle[i],'NB']
        df.append(da1)
        
        lins = []
        inp2 = os.path.join(dirname,file2[i])
        with open(inp2,'r') as f:
            lines = f.read().splitlines()
            for line in lines:
                nm = line.split(' ')
                ax = list(filter(None,nm))
                lins.append(ax)

            da = pd.DataFrame(lins)  
            da1 = da.iloc[13:,:]
            da2 = da1.drop(da.tail(3).index)
            da2.index = np.arange(da2.shape[0])
            dal1 = da2.iloc[:,[2,11]]
            dal1.columns = ['Peptide',alle[i]]
            df2.append(dal1)
            
    dax = df[0]
    dalx = df2[0]
    if len(alle)>=2 :
        for i in range(1,len(alle)):
            dax = pd.merge(dax,df[i],on='Peptide')
            dalx = pd.merge(dalx,df2[i],on='Peptide')
    return dax,dalx
    
def calculate(dax,dalx,alle,dap):
    
    dalx[alle] = dalx[alle].apply(pd.to_numeric)
    dalx1 = dalx.loc[:,alle]

    x = dalx1.T.apply(lambda x: x.nsmallest(2)).T
    dalx['annotation_score'] = round(x.max(axis=1)/x.min(axis=1)).astype(int)

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

    dax4 = pd.merge(dax3,dalx,how='outer',on=['Peptide'])
    dax4.columns = dax4.columns.str.replace('Peptide','Peptide_Sequence')
    
    dan = pd.merge(dap,dax4,how='outer',on=['Peptide_Sequence'])
    dan['top_allele']= dan.loc[:,'top_allele'].fillna('unclassified')
    
    return dan

def plotheatmap(dalx,dirname):
    
    dalx[alle] = dalx[alle].apply(pd.to_numeric)
    dalx_temp = dalx[alle].copy()
    dalx_temp[dalx_temp > 1000] = 1000
    rsp_df = dalx_temp
    rsp_df.sort_values(alle, inplace=True)
    rsp_df.index = np.arange(rsp_df.shape[0])
    rsp_df

    sns.set(font_scale=1.5)
#plt.rc('font',family='Times New Roman',size=18)
    f, ax = plt.subplots(figsize = (15,12))
    sns.heatmap(rsp_df, fmt="g", cmap='coolwarm')
    ax.invert_yaxis()
#sns.heatmap(wida2,vmax=500,vmin=40,center=True,fmt='d',annot=True,ax=ax,cmap='coolwarm',annot_kws={'size':18})
    plt.setp(ax.get_xticklabels(), rotation=30,fontsize=15, horizontalalignment='right')

    lab = []
    if int(rsp_df.shape[0]/10000)!=0 :
        num_tick = int(rsp_df.shape[0]/10000)+1
        for i in range(num_tick+1):
            lab.append(10000*i)
    elif int(rsp_df.shape[0]/1000)!=0 :
        num_tick = int(rsp_df.shape[0]/1000)+1
        for i in range(num_tick+1):
            lab.append(1000*i)
    elif int(rsp_df.shape[0]/100)!=0 :
        num_tick = int(rsp_df.shape[0]/100)+1
        for i in range(num_tick+1):
            lab.append(100*i)
    elif int(rsp_df.shape[0]/10)!=0 :
        num_tick = int(rsp_df.shape[0]/10)+1
        for i in range(num_tick+1):
            lab.append(10*i)
        
    ax.set_yticks(lab)
    ax.set_yticklabels(lab)

    ax.set_ylabel('Peptide Counts',fontsize=20)
    plt.tight_layout()
    out1 = 'annotation_heatmap.svg'
    plt.savefig(out1,dpi=800)
    
    return 'ok'
#plt.show()



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
            if 'xls' in i:
                file.append(i)
        file1 = sorted(file)
        for i in file1:
            dirStr, ext = os.path.splitext(i)
            
            a1 = re.sub('xls','log',i)
            file2.append(a1)
            if 'D' in dirStr:
                a = dirStr
            else:
                a=re.sub('_',':',dirStr)
            alle.append(a)
             
        dapx = pd.read_csv(infile)

        if 'D' in alle[0]:
            dax,dalx = readc2(file1,alle,file2)
            cla = 'classII'
            dap = dapx.loc[(dapx['Peptide_Sequence'].str.len()>=9) & (dapx['Peptide_Sequence'].str.len()<=25)]
        else:
            dax,dalx = readc1(file1,alle,file2)
            cla = 'classI'
            dap = dapx.loc[(dapx['Peptide_Sequence'].str.len()>=8) & (dapx['Peptide_Sequence'].str.len()<=14)]

        plotheatmap(dalx,dirname)
        
        danx = calculate(dax,dalx,alle,dap)
        dan = danx.copy()
        dan['length'] = dan['Peptide_Sequence'].str.len()
        out3 = cla + '_peptide_complete.csv'
        dan.to_csv(out3,index=False)

        dan_temp = dan[alle]
        out2 = 'alle_annotation.csv'
        dan_temp.to_csv(out2,index=False)
