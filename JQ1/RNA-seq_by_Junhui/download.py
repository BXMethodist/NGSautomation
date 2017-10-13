# import os
#
# f = open('Samples.txt', 'r')
# results = []
# for line in f:
#     line = line.strip()
#     if line.endswith('.gz'):
#         link = 'ftp://128.120.88.242/data_release/C202SC17071324/raw_data/'+ line
#         results.append(line)
#         os.system('wget --user P202SC17071641-01_20170823_VSY0Hs --password ZcUjLj '+link)
#
# f = open('Samples.txt', 'w')
# for r in results:
#     f.write(line+'\n')
# f.close()

import pandas as pd, numpy as np
# import matplotlib.pyplot as plt
#
# df = pd.read_csv('MQ_CT_FC2.csv', index_col=0)
# df['M'] = pd.to_numeric(df['M'], errors='coerce')
# df['tissue'] = pd.to_numeric(df['tissue'], errors='coerce')
#
# df = df.fillna(-np.inf)
#
# df.plot()
# # -0.042203
# df.plot(x='M', y='tissue', kind='scatter')
# plt.show()

df = pd.read_csv('./cuffdiff/gene_exp.diff',sep='\t')
print df['sample_1'].unique()
print df['sample_2'].unique()

for group1 in df['sample_1'].unique():
    for group2 in df['sample_2'].unique():
        cur_df = df[(df['sample_1']==group1) & (df['sample_2']==group2)]
        cur_df.to_csv('./cuffdiff/'+group1+'_'+group2+'_gene_exp.diff.tsv', sep='\t')