#!/usr/bin/env python3
#ptm_heatmap.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

print ('usage: please input argv "argv1:modify.csv')

data1 = pd.read_csv(sys.argv[1],header=None,sep='\s+')
daname = pd.read_csv(sys.argv[2],header=None)

data1.columns = ['ions','peptide','charge','aminoa','ptm','location']
data1['num'] = 1
data = data1.copy()
df1 = data.drop(data[data['ptm'].str.contains('Homoserine',na=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('USM',na=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('Met->Glu',na=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('Cys->Thr',na=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('2H(4)',na=False,regex=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('hydroxymethyl',na=False,regex=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('Val',na=False,regex=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('Acetyl13C(6)-Silac-label',na=False,regex=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('Asp->Tyr',na=False,regex=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('Propyl:2H(6)',na=False,regex=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('MDA-adduct+62',na=False,regex=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('Gly->His',na=False,regex=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('TMT6plex',na=False,regex=False)].index)
df1 = df1.drop(df1[df1['ptm'].str.contains('Delta:S(-1)Se(1)',na=False,regex=False)].index)

old_name = list(daname.iloc[:,1])
new_name = list(daname.iloc[:,0])

df_new = df1['ptm'].replace(old_name,new_name)
data4 = df1.copy()
data4['ptm1'] = df_new

data4.to_csv('modify_info.csv',index=False)

da = data4.loc[:,['aminoa','ptm1','num']]
wida = da.pivot_table(index='aminoa', columns='ptm1',values='num',aggfunc=[np.sum],fill_value=0)
wida1 = pd.DataFrame(wida.values.T,index=wida.columns, columns=wida.index)

ind=[]
for i in range(len(wida1.index)):
    ind.append(wida1.index[i][1])
wida2 = pd.DataFrame(wida1.values,index=ind,columns=wida1.columns)

wida2['total'] = wida2.apply(lambda x: x.sum(),axis=1)
wida2.sort_values(by='total',ascending=False,inplace=True)
wida2.loc['sum'] = wida2.apply(lambda x: x.sum())
#wida3=wida2.sort_values(by='sum',axis=1,ascending=True)

sns.set(font_scale=1.5)
plt.rc('font',family='Times New Roman',size=18)
wida2.iloc[0:-1,:-1].plot.barh(stacked=True,figsize=(15,15))
plt.title('Modifications')
plt.xlabel('Number')
plt.ylabel('PTMs')
plt.xscale('symlog')
plt.tight_layout()
plt.savefig('barplot_ptm.svg',dpi=600)
plt.savefig('barplot_ptm.png',dpi=600)
#plt.show()

sns.set(font_scale=1.5)
plt.rc('font',family='Times New Roman',size=18)
f, ax = plt.subplots(figsize = (20,16))
sns.heatmap(wida2,vmax=30,vmin=-30,fmt='d',annot=True,cmap='RdBu_r',ax=ax,annot_kws={'size':18})
plt.setp(ax.get_xticklabels(), rotation=-30,fontsize=18, horizontalalignment='right')
plt.setp(ax.get_yticklabels(),fontsize=20)
#ax.set_title('Colormap',fontsize=28)
ax.set_xlabel('Amino Acids',fontsize=25)
ax.set_ylabel('PTMs',fontsize=25)
plt.tight_layout()
plt.savefig('heatmap_ptm.svg',dpi=600)
plt.savefig('heatmap_ptm.png',dpi=600)
#plt.show()


