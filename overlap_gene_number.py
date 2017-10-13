import pandas as pd, os, numpy as np

# gtf = pd.read_csv('./JQ1/mm9.20150218.knownGene.xls', sep='\t')
#
# tables = ['./JQ1/unique/' + x for x in os.listdir('./JQ1/unique/') if x.endswith('_gtf.xls')]
# gene_table = ['./JQ1/unique/' + x for x in os.listdir('./JQ1/unique/') if x.endswith('.txt') and x.startswith('M')]
#
# t = './JQ1/unique/down_SRR3929171_total_sig.txtgene_gtf.xls'
# k = './JQ1/unique/M0_down_list.txt'
# print t, k
# df = pd.read_csv(t, sep='\t')
# exp_df = pd.read_csv(k, sep='\t', header=None)
# super_genes = set(df['mm9.kgXref.geneSymbol'].unique())
# exp_genes = set(exp_df.ix[:, 0].unique())
# print len(super_genes.intersection(exp_genes))
# print len(super_genes), len(exp_genes)
# print len(super_genes)*len(exp_genes)*1.0/29955
#
#
# t = './JQ1/unique/up_SRR3929171_total_sig.txtgene_gtf.xls'
# k = './JQ1/unique/M0_up_list.txt'
# print t, k
# df = pd.read_csv(t, sep='\t')
# exp_df = pd.read_csv(k, sep='\t', header=None)
# super_genes = set(df['mm9.kgXref.geneSymbol'].unique())
# exp_genes = set(exp_df.ix[:, 0].unique())
# print len(super_genes.intersection(exp_genes))
# print len(super_genes), len(exp_genes)
# print len(super_genes)*len(exp_genes)*1.0/29955
#
# t = './JQ1/unique/down_SRR3929176_total_sig.txtgene_gtf.xls'
# k = './JQ1/unique/M1_down_list.txt'
# print t, k
# df = pd.read_csv(t, sep='\t')
# exp_df = pd.read_csv(k, sep='\t', header=None)
# super_genes = set(df['mm9.kgXref.geneSymbol'].unique())
# exp_genes = set(exp_df.ix[:, 0].unique())
# print len(super_genes.intersection(exp_genes))
# print len(super_genes), len(exp_genes)
# print len(super_genes)*len(exp_genes)*1.0/29955
#
# t = './JQ1/unique/up_SRR3929176_total_sig.txtgene_gtf.xls'
# k = './JQ1/unique/M1_up_list.txt'
# print t, k
# df = pd.read_csv(t, sep='\t')
# exp_df = pd.read_csv(k, sep='\t', header=None)
# super_genes = set(df['mm9.kgXref.geneSymbol'].unique())
# exp_genes = set(exp_df.ix[:, 0].unique())
# print super_genes.intersection(exp_genes)
# print len(super_genes), len(exp_genes)
# print len(super_genes)*len(exp_genes)*1.0/29955
#
#
#
# t = './JQ1/unique/down_SRR3929175_total_sig.txtgene_gtf.xls'
# k = './JQ1/unique/M2_down_list.txt'
# print t, k
# df = pd.read_csv(t, sep='\t')
# exp_df = pd.read_csv(k, sep='\t', header=None)
# super_genes = set(df['mm9.kgXref.geneSymbol'].unique())
# exp_genes = set(exp_df.ix[:, 0].unique())
# print len(super_genes.intersection(exp_genes))
# print len(super_genes), len(exp_genes)
# print len(super_genes)*len(exp_genes)*1.0/29955
#
# t = './JQ1/unique/up_SRR3929175_total_sig.txtgene_gtf.xls'
# k = './JQ1/unique/M2_up_list.txt'
# print t, k
# df = pd.read_csv(t, sep='\t')
# exp_df = pd.read_csv(k, sep='\t', header=None)
# super_genes = set(df['mm9.kgXref.geneSymbol'].unique())
# exp_genes = set(exp_df.ix[:, 0].unique())
# print super_genes.intersection(exp_genes)
# print len(super_genes), len(exp_genes)
# print len(super_genes)*len(exp_genes)*1.0/29955
#
# print len(gtf['mm9.kgXref.geneSymbol'].unique())


CT_df = pd.read_csv('./JQ1/RNA-seq_by_Junhui/cuffdiff/C_T_significant.txt', sep='\t')
MQ_df = pd.read_csv('./JQ1/RNA-seq_by_Junhui/cuffdiff_MQ/M_Q_avg.below.remove_2fold.txt', sep='\t')

down_CT = CT_df[CT_df['log2(fold_change)']=='#NAME?']['gene_id'].tolist()
CT_df = CT_df[CT_df['log2(fold_change)']!='#NAME?']
CT_df['log2(fold_change)'] = pd.to_numeric(CT_df['log2(fold_change)'])
down_CT += CT_df[CT_df['log2(fold_change)']<0]['gene_id'].tolist()
up_CT = CT_df[CT_df['log2(fold_change)']>0]['gene_id'].tolist()

down_MQ = MQ_df[MQ_df['log2(fold_change)']=='#NAME?']['gene_id'].tolist()
MQ_df = MQ_df[MQ_df['log2(fold_change)']!='#NAME?']
MQ_df['log2(fold_change)'] = pd.to_numeric(MQ_df['log2(fold_change)'])
down_MQ += MQ_df[MQ_df['log2(fold_change)']<0]['gene_id'].tolist()
up_MQ = MQ_df[MQ_df['log2(fold_change)']>0]['gene_id'].tolist()

# print len(down_CT), len(up_CT), len(down_MQ), len(up_MQ)
down_CT = set(down_CT)
up_CT = set(up_CT)
down_MQ = set(down_MQ)
up_MQ = set(up_MQ)

print len(down_CT.intersection(down_MQ))
print len(down_CT)-len(down_CT.intersection(down_MQ)), len(down_MQ)-len(down_CT.intersection(down_MQ))
print len(down_CT)*len(down_MQ)*1.0/29955

print len(up_CT.intersection(up_MQ))
print len(up_CT)-len(up_CT.intersection(up_MQ)), len(up_MQ)-len(up_CT.intersection(up_MQ))
print len(up_CT)*len(up_MQ)*1.0/29955

print down_CT.intersection(down_MQ)

print up_CT.intersection(up_MQ)