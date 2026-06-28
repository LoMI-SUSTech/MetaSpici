#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import sys


# In[2]:


fn1 = sys.argv[1]
#fn1="bin2855_sel.bed"


# In[3]:


raw_bed_file = pd.read_csv(fn1, sep='\t', header=None)
bed_file=raw_bed_file[[0,1,2,3]]
bed_file.rename(columns={0:"bin_orf",1:"start_loc",2:"end_loc",3:"reads"}, inplace=True)

bed_file['read_name']=bed_file["reads"].str.split(":",n=0,expand=True)[0]


# In[6]:


merge_bed = pd.merge(bed_file, bed_file, on="read_name",how="inner")

merge_bed['type1'] = np.where(merge_bed.eval('reads_x==reads_y'),'True','False')
merge_bed['type2'] = np.where(merge_bed.eval('bin_orf_x==bin_orf_y'),'True','False')

new_merge_bed = merge_bed[merge_bed['type1']== "False"]

new_merge_bed.drop_duplicates(subset="read_name",inplace=True)
new_merge_bed['sample'] = str(fn1).split("_")[0]
new_merge_bed[new_merge_bed['type2']== 'False'].to_csv(fn1+".ciro_diff.txt",sep='\t',index=0)
new_merge_bed[new_merge_bed['type2']== 'True'].to_csv(fn1+".ciro_same.txt",sep='\t',index=0)

