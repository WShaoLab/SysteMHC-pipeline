#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import re

#print ('usage: please input argv "argv1:peptidelist", "argv2:ionslist","argv3:MHCpan.xls","argv4:outname"')

data = pd.read_csv(sys.argv[2])
datare = pd.read_csv(sys.argv[1])
alle = sys.argv[3]

data.columns = data.columns.str.replace('search_hit','Peptide_Sequence')

dax = data[data['top_allele']==alle]
pep = dax.loc[:,'Peptide_Sequence']
peplist = pep.tolist()
dataions = datare[datare['Peptide_Sequence'].isin(peplist)]
bindions = dataions.loc[:,'Ions']

aname = re.sub(':','_',alle)

outn1 = aname + '_bindions.csv'
outn2 = aname + '_bindpeps.csv'
    
bindions.to_csv(outn1,header=None,index=False)
pep.to_csv(outn2,header=None,index=False)

