import os, pandas as pd

def cuffdiff_table_selector(table, name, **kwargs):
    ## **kwargs pair of column name and cutoffs
    ##TO DO: create table selecter method for heatmaps and downstream analysis
    df = pd.read_csv(table, sep="\t", index_col="gene_id")
    output_names = [name]
    result_df = None

    for key, value in kwargs.items():
        if key != "gene_id":
            output_names.append(key)
            output_names.append(str(value))
        if result_df is None:
            if key == "gene_id":
                print value
                result_df = df.ix[value, :]
                print result_df
            if key == "q_value":
                result_df = df[df['q_value'] < value]
            if key == "p_value":
                result_df = df[df['p_value'] < value]
            if key == "log2(fold_change)":
                if value > 0:
                    result_df = df[df['log2(fold_change)'] > value]
                if value < 0:
                    result_df = df[df['log2(fold_change)'] < value]
            if key == "sample_1":
                result_df = df[df['sample_1'] == value]
            if key == "sample_2":
                result_df = df[df['sample_2'] == value]
        else:
            if key == "gene_id":
                result_df = result_df.ix[value, :]
            if key == "q_value":
                result_df = result_df[result_df['q_value'] < value]
            if key == "p_value":
                result_df = result_df[result_df['p_value'] < value]
            if key == "log2(fold_change)":
                if value > 0:
                    result_df = result_df[result_df['log2(fold_change)'] > value]
                if value < 0:
                    result_df = result_df[result_df['log2(fold_change)'] < value]
            if key == "sample_1":
                result_df = result_df[result_df['sample_1'] == value]
            if key == "sample_2":
                result_df = result_df[result_df['sample_2'] == value]

    output_filename = "_".join(output_names) + ".tsv"

    del result_df['test_id']
    del result_df['gene']
    del result_df['locus']
    del result_df['test_stat']

    result_df.to_csv(output_filename, sep="\t")

    return result_df, "_".join(output_names)

f = open('genelist_gianni', "r")

gene = []

for line in f.readlines():
    gene.append(line.strip())

f.close()

char = {'gene_id': gene}

df1 = cuffdiff_table_selector('gene_exp.diff','Gianni_interested_genes_values.tsv', **char)[0]
