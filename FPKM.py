import pandas as pd, numpy as np

def scale(y, c=True, sc=True):
    x = y.copy()

    if c:
        x -= x.mean()
    if sc and c:
        x /= x.std()
    elif sc:
        x /= np.sqrt(x.pow(2).sum().div(x.count() - 1))
    return x

def get_fpkm(significant_gene_table_path, cuff_diff_fpkm_table_path, groupnameA, groupnameB):
    FPKM_df = pd.read_csv(cuff_diff_fpkm_table_path, sep='\t')
    all_replicates = sorted(FPKM_df['replicate'].unique())
    expression = {}
    for i in range(FPKM_df.shape[0]):
        gene = FPKM_df.ix[i, 'tracking_id']
        condition = FPKM_df.ix[i, 'condition']
        rep = FPKM_df.ix[i, 'replicate']
        expr = FPKM_df.ix[i, 'FPKM']
        expression[(gene, condition, rep)] = expr
    print 'load complete'

    gene_df = pd.read_csv(significant_gene_table_path, sep='\t')

    genes = gene_df['gene_id'].unique()

    final = []
    columns = []

    for sample in [groupnameA, groupnameB]:
        for rep in all_replicates:
            columns.append(sample+'_'+str(rep))

    for gene in genes:
        cur_final = []
        for sample in [groupnameA, groupnameB]:
            for rep in all_replicates:
                if (gene, sample, rep) in expression:
                    cur_final.append(expression[(gene, sample, rep)])
        # cur_final = scale(np.asarray(cur_final))
        cur_final = [gene] + list(cur_final)
        final.append(cur_final)

    final_df = pd.DataFrame(final)
    final_df.columns = ['NAME']+columns
    final_df['NAME']=final_df['NAME'].str.upper()
    final_df.to_csv(groupnameA+groupnameB+'_fpkm.tsv', index=None, sep='\t')
    return

# get_fpkm('./JQ1/RNA-seq_by_Junhui/cuffdiff/gene_exp.diff','./JQ1/RNA-seq_by_Junhui/cuffdiff/genes.read_group_tracking', 'C','T')
# get_fpkm('./JQ1/RNA-seq_by_Junhui/cuffdiff/gene_exp.diff','./JQ1/RNA-seq_by_Junhui/cuffdiff/genes.read_group_tracking', 'Cn','Tn')
#
# get_fpkm('./JQ1/RNA-seq_by_Junhui/cuffdiff_MQ/gene_exp.diff','./JQ1/RNA-seq_by_Junhui/cuffdiff_MQ/genes.read_group_tracking', 'M','Q')
# get_fpkm('./JQ1/RNA-seq_by_Junhui/cuffdiff/Cn_Tn_significant.txt','./JQ1/RNA-seq_by_Junhui/cuffdiff/genes.read_group_tracking', 'Cn','Tn')
get_fpkm('./PAF/novogene/cuffdiff/gene_exp.diff','./PAF/novogene/cuffdiff/genes.read_group_tracking', 'WT','PAF')