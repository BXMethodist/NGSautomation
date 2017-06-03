import pandas as pd

df = pd.read_csv('PIC0_30.csv', index_col=0)

gene_ids = []

gene_df = pd.read_csv('mouse_histone.txt', sep='\t', index_col=0)

for i in range(df.shape[0]):
    gene_id = df.ix[i, 'gene_id']
    gene_id = gene_id.lower().strip()
    gene_ids.append(gene_id)

df['gene_id'] = gene_ids
df = df.set_index(['gene_id'])
genes = [str(x).lower().strip() for x in gene_df.index]
print genes
sub_df = df.ix[genes, :]

print sub_df.shape

sub_df.to_csv('PIC0_30_histone.csv')