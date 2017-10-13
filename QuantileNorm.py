import pandas as pd

df = pd.read_csv('./PAF/genes_read_group_FPKM_GSEA.txt', sep='\t', dtype={"NAME":str})

df = df.set_index(['NAME'])

rownames = set()
dups = set()
for i in df.index:
    if i not in rownames:
        rownames.add(i)
    else:
        dups.add(i)
df = df[~df.index.isin(dups)]

rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()

result_df = df.rank(method='min').stack().astype(int).map(rank_mean).unstack()

print result_df

result_df.to_csv('./PAF/novogene/genes_read_group_FPKM_GSEA_QT.txt')