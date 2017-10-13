import pandas as pd, numpy as np
from collections import defaultdict

# def scale(y, c=True, sc=True):
#     x = y.copy()
#
#     if c:
#         x -= x.mean()
#     if sc and c:
#         x /= x.std()
#     elif sc:
#         x /= np.sqrt(x.pow(2).sum().div(x.count() - 1))
#     return x
#
# df = pd.read_csv('genes.read_group_tracking', sep='\t')
#
# # print df
# results = defaultdict(float)
# for i in range(df.shape[0]):
#     gene = df.ix[i, 'tracking_id']
#     condition = df.ix[i, 'condition']
#     replicate = df.ix[i, 'replicate']
#     FPKM = df.ix[i, 'FPKM']
#     results[(gene, condition, replicate)] = FPKM
#
# samples = [(condition, replicate) for condition in df['condition'].unique()
#                                   for replicate in df['replicate'].unique()]
#
# ###
# indexes = list(pd.read_csv('./gene_exp_diff_significant.txt', sep='\t')['gene_id'].unique())
# final = []
# for gene in df['tracking_id'].tolist():
#     ####
#     if gene in indexes:
#         cur_final = [results[tuple([gene]+list(sample))] for sample in samples]
#         cur_final = scale(np.asarray(cur_final))
#         final.append([gene]+list(cur_final))
#
# names = ['_'.join([str(i) for i in s]) for s in samples]
# final_df = pd.DataFrame(final)
# final_df.columns = ['NAME']+names
#
# final_df = final_df.drop_duplicates()

# final_df['NAME'] = final_df['NAME'].str.upper()
# print final_df['NAME']

# final_df.to_csv('genes_read_group_FPKM_GSEA.txt', sep='\t', index=None)

# indexes = list(pd.read_csv('./gene_exp_diff_significant.txt', sep='\t')['gene_id'].unique())


# final_df = final_df[final_df['NAME'].isin(indexes)]
# final_df.to_csv('genes_read_group_FPKM_significant.txt', sep='\t', index=None)

df = pd.read_csv('log2Rank.txt', sep='\t', index_col=0)
# print df
df.index = df.index.str.upper()

df.to_csv('RankbyQ.rnk', sep='\t')