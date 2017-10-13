## Three way overlap 6 groups of genes

import pandas as pd
import numpy as np, os
from copy import deepcopy

def scale(y, c=True, sc=True):
    x = y.copy()

    if c:
        x -= x.mean()
    if sc and c:
        x /= x.std()
    elif sc:
        x /= np.sqrt(x.pow(2).sum().div(x.count() - 1))
    return x

def to_list(gene_list, file_name):
    f = open(file_name, 'w')
    for g in gene_list:
        f.write(g+'\n')
    f.close()

# gtf = pd.read_csv('./JQ1/mm9.20150218.knownGene.xls', sep='\t')
#
# df = pd.read_csv('./JQ1/gene_exp.diff_significant.tsv', sep='\t',
#                  dtype={'log2(fold_change)': np.float64} )
#
# pairs = [['M0_UT', 'M1_IFN', 'M2_IL4'],
#          ['M1_IFN', 'M0_UT', 'M2_IL4'],
#          ['M2_IL4', 'M0_UT','M1_IFN'],
#          ]
#
# results = []
#
#
# s1, s2, s3 = 'M0_UT', 'M1_IFN', 'M2_IL4'
#
# ### specific up or down in M0_UT
# df1 = df[(df['sample_1'] == s1) & (df['sample_2'] == s2)]
# df2 = df[(df['sample_1'] == s1) & (df['sample_2'] == s3)]
# #
# up_df1 = df1[df1['log2(fold_change)'] > 0]
# up_df2 = df2[df2['log2(fold_change)'] > 0]
#
# cur_results = set()
#
# for gene in up_df1['gene_id'].values:
#     if gene in up_df2['gene_id'].values:
#         cur_results.add(gene)
#
# results.append(cur_results)
#
# to_list(cur_results, './JQ1/M0_down_list.txt')
# print len(cur_results)
#
# for g in cur_results:
#     if g not in gtf['mm9.kgXref.geneSymbol'].unique():
#         print g
#
# cur_gtf = gtf[gtf['mm9.kgXref.geneSymbol'].isin(cur_results)]
# # cur_gtf = cur_gtf[['mm9.knownGene.chrom',	'mm9.knownGene.txStart', 'mm9.knownGene.txEnd']].drop_duplicates()
# cur_gtf.to_csv('./JQ1/M0_down_list.tsv', sep='\t',index=None, header=None)
#
#
# down_df1 = df1[df1['log2(fold_change)'] < 0]
# down_df2 = df2[df2['log2(fold_change)'] < 0]
#
# cur_results = set()
#
# # print down_df1.shape, down_df2.shape
#
# for gene in down_df1['gene_id'].values:
#     if gene in down_df2['gene_id'].values:
#         cur_results.add(gene)
#
# print len(cur_results.intersection(BRD4_genes))
# print len(cur_results)
#
# to_list(cur_results, './JQ1/M0_up_list.txt')
# results.append(cur_results)
#
# for g in cur_results:
#     if g not in gtf['mm9.kgXref.geneSymbol'].unique():
#         print g
#
# cur_gtf = gtf[gtf['mm9.kgXref.geneSymbol'].isin(cur_results)]
# # cur_gtf = cur_gtf[['mm9.knownGene.chrom',	'mm9.knownGene.txStart', 'mm9.knownGene.txEnd']].drop_duplicates()
# cur_gtf.to_csv('./JQ1/M0_up_list.tsv', sep='\t',index=None, header=None)
#
# # print [len(x) for x in results]
#
# ### specific up or down in M1_IFN
# df1 = df[(df['sample_1'] == s1) & (df['sample_2'] == s2)]
# df2 = df[(df['sample_1'] == s2) & (df['sample_2'] == s3)]
# #
# up_df1 = df1[df1['log2(fold_change)'] < 0]
# up_df2 = df2[df2['log2(fold_change)'] > 0]
#
# cur_results = set()
#
# for gene in up_df1['gene_id'].values:
#     if gene in up_df2['gene_id'].values:
#         cur_results.add(gene)
#
# results.append(cur_results)
# to_list(cur_results, './JQ1/M1_down_list.txt')
# for g in cur_results:
#     if g not in gtf['mm9.kgXref.geneSymbol'].unique():
#         print g
#
# cur_gtf = gtf[gtf['mm9.kgXref.geneSymbol'].isin(cur_results)]
# # cur_gtf = cur_gtf[['mm9.knownGene.chrom',	'mm9.knownGene.txStart', 'mm9.knownGene.txEnd']].drop_duplicates()
# cur_gtf.to_csv('./JQ1/M1_down_list.tsv', sep='\t',index=None, header=None)
#
#
# down_df1 = df1[df1['log2(fold_change)'] > 0]
# down_df2 = df2[df2['log2(fold_change)'] < 0]
#
# cur_results = set()
#
# # print down_df1.shape, down_df2.shape
#
# for gene in down_df1['gene_id'].values:
#     if gene in down_df2['gene_id'].values:
#         cur_results.add(gene)
#
# print len(cur_results.intersection(BRD4_genes))
# print len(cur_results)
# results.append(cur_results)
# to_list(cur_results, './JQ1/M1_up_list.txt')
#
# for g in cur_results:
#     if g not in gtf['mm9.kgXref.geneSymbol'].unique():
#         print g
#
# cur_gtf = gtf[gtf['mm9.kgXref.geneSymbol'].isin(cur_results)]
# # cur_gtf = cur_gtf[['mm9.knownGene.chrom',	'mm9.knownGene.txStart', 'mm9.knownGene.txEnd']].drop_duplicates()
# cur_gtf.to_csv('./JQ1/M1_up_list.tsv', sep='\t',index=None, header=None)
#
# # print [len(x) for x in results]
#
# ### specific up or down in M2_IL4
# df1 = df[(df['sample_1'] == s1) & (df['sample_2'] == s3)]
# df2 = df[(df['sample_1'] == s2) & (df['sample_2'] == s3)]
# #
# up_df1 = df1[df1['log2(fold_change)'] < 0]
# up_df2 = df2[df2['log2(fold_change)'] < 0]
#
# cur_results = set()
#
# for gene in up_df1['gene_id'].values:
#     if gene in up_df2['gene_id'].values:
#         cur_results.add(gene)
#
# to_list(cur_results, './JQ1/M2_down_list.txt')
# results.append(cur_results)
#
# for g in cur_results:
#     if g not in gtf['mm9.kgXref.geneSymbol'].unique():
#         print g
#
# cur_gtf = gtf[gtf['mm9.kgXref.geneSymbol'].isin(cur_results)]
# # cur_gtf = cur_gtf[['mm9.knownGene.chrom',	'mm9.knownGene.txStart', 'mm9.knownGene.txEnd']].drop_duplicates()
# cur_gtf.to_csv('./JQ1/M2_down_list.tsv', sep='\t',index=None, header=None)
#
#
# down_df1 = df1[df1['log2(fold_change)'] > 0]
# down_df2 = df2[df2['log2(fold_change)'] > 0]
#
# cur_results = set()
#
# # print down_df1.shape, down_df2.shape
#
# for gene in down_df1['gene_id'].values:
#     if gene in down_df2['gene_id'].values:
#         cur_results.add(gene)
# print len(cur_results.intersection(BRD4_genes))
# print len(cur_results)
# results.append(cur_results)
# to_list(cur_results, './JQ1/M2_up_list.txt')
#
# for g in cur_results:
#     if g not in gtf['mm9.kgXref.geneSymbol'].unique():
#         print g
#
# cur_gtf = gtf[gtf['mm9.kgXref.geneSymbol'].isin(cur_results)]
# # cur_gtf = cur_gtf[['mm9.knownGene.chrom',	'mm9.knownGene.txStart', 'mm9.knownGene.txEnd']].drop_duplicates()
# cur_gtf.to_csv('./JQ1/M2_up_list.tsv', sep='\t',index=None, header=None)


# print [len(x) for x in results]

# print scale(np.asarray([1.0,2.0,3.0,4.0]))

# FPKM_df = pd.read_csv('./JQ1/genes.read_group_tracking', sep='\t')
# all_samples = sorted(FPKM_df['condition'].unique())
# all_replicates = sorted(FPKM_df['replicate'].unique())
#
# print all_samples, all_replicates
#
# expression = {}
# for i in range(FPKM_df.shape[0]):
#     gene = FPKM_df.ix[i, 'tracking_id']
#     condition = FPKM_df.ix[i, 'condition']
#     rep = FPKM_df.ix[i, 'replicate']
#     expr = FPKM_df.ix[i, 'FPKM']
#     expression[(gene, condition, rep)] = expr
# print 'load complete'
# final = []
# # columns = []
# for group in results:
#     print len(group)
#     for gene in group:
#         cur_final = []
#         for sample in all_samples:
#             for rep in all_replicates:
#                 # print cur_FPKM_df.columns
#                 if (gene, sample, rep) in expression:
#                     cur_final.append(expression[(gene, sample, rep)])
#         cur_final = scale(np.asarray(cur_final))
#         cur_final = [gene] + list(cur_final)
#         if len(cur_final) == 10:
#             final.append(cur_final)
#     final.append([' ']*10)
# final_df = pd.DataFrame(final)
#
# final_df.to_csv('./JQ1/gene_unique_exp_diff.tsv', sep='\t', index=None)


### for super enhance unique

results = []

tables = ['./JQ1/unique/'+x for x in os.listdir('./JQ1/unique') if x.endswith('_gtf.xls')]
print tables
for t in tables:
    t_df = pd.read_csv(t, sep='\t')
    results.append(set(t_df['mm9.kgXref.geneSymbol'].unique()))
    print len(set(t_df['mm9.kgXref.geneSymbol'].unique()))



FPKM_df = pd.read_csv('./JQ1/genes.read_group_tracking', sep='\t')
all_samples = sorted(FPKM_df['condition'].unique())
all_replicates = sorted(FPKM_df['replicate'].unique())

print all_samples, all_replicates

expression = {}
for i in range(FPKM_df.shape[0]):
    gene = FPKM_df.ix[i, 'tracking_id']
    condition = FPKM_df.ix[i, 'condition']
    rep = FPKM_df.ix[i, 'replicate']
    expr = FPKM_df.ix[i, 'FPKM']
    expression[(gene, condition, rep)] = expr
print 'load complete'
final = []
# columns = []
for group in results:
    print len(group)
    for gene in group:
        cur_final = []
        for sample in all_samples:
            for rep in all_replicates:
                # print cur_FPKM_df.columns
                if (gene, sample, rep) in expression:
                    cur_final.append(expression[(gene, sample, rep)])
        cur_final = scale(np.asarray(cur_final))
        cur_final = [gene] + list(cur_final)
        if len(cur_final) == 10:
            final.append(cur_final)
    final.append([' ']*10)
final_df = pd.DataFrame(final)

# final_df.to_csv('./JQ1/superenhancer_gene_unique_exp_diff.tsv', sep='\t', index=None)