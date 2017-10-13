import pandas as pd
from collections import defaultdict

df = pd.read_csv('genes.read_group_tracking', sep='\t')

# print df
results = defaultdict(float)
for i in range(df.shape[0]):
    gene = df.ix[i, 'tracking_id']
    condition = df.ix[i, 'condition']
    replicate = df.ix[i, 'replicate']
    FPKM = df.ix[i, 'FPKM']
    results[(gene, condition, replicate)] = FPKM

samples = [(condition, replicate) for condition in df['condition'].unique()
                                  for replicate in df['replicate'].unique()]

final = []
for gene in df['tracking_id'].tolist():
    final.append([gene]+[results[tuple([gene]+list(sample))] for sample in samples])

names = ['_'.join([str(i) for i in s]) for s in samples]
final_df = pd.DataFrame(final)
final_df.columns = ['NAME']+names

final_df = final_df.drop_duplicates()

final_df['NAME'] = final_df['NAME'].str.upper()
# print final_df['NAME']

final_df.to_csv('genes_read_group_FPKM_GSEA.txt', sep='\t', index=None)
