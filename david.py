import os, pandas as pd, numpy as np

def convert_david(path, type):
    df = pd.read_csv(path, sep='\t')
    df = df[['Term', 'Count', '%', 'Genes', 'Fold Enrichment', 'Benjamini']]
    df['-log10 enrich P'] = -np.log10(df['Benjamini'])
    if type == 'KEGG':
        new_terms = [term[1] for term in df['Term'].str.split(':').values]
        df['Term'] = new_terms
    elif type == 'GO':
        new_terms = [term[1] for term in df['Term'].str.split('~').values]
        df['Term'] = new_terms

    df.to_csv(path, sep='\t', index=None)




path ='./JQ1/RNA-seq_by_Junhui/David/M_Q_up_in_Q_KEGG.txt'
type = 'KEGG'

convert_david(path, type)

path ='./JQ1/RNA-seq_by_Junhui/David/M_Q_up_in_Q_GO.BP.txt'
type = 'GO'

convert_david(path, type)
