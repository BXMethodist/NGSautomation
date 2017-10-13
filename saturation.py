import pandas as pd
from bisect import bisect_left

def df_to_index_danpos(df, bin=1000):
    results = {}
    results_indexes = {}
    f = open(df, 'r')

    for line in f.readlines()[1:]:
        line = line.split()
        if float(line[5]) == 0:
            continue
        t = (line[0],
             int(line[1]),
             int(line[2]),
             int(line[4]),
             float(line[5]),
             float(line[6]),
             )

        if t[0] not in results:
            results[t[0]] = {}

        for i in range(t[1]/bin, t[2]/bin+1):
            if i in results[t[0]]:
                results[t[0]][i].add(t)
            else:
                results[t[0]][i] = set()
                results[t[0]][i].add(t)
    f.close()

    for key in results.keys():
        cur_indexes = results[key].keys()
        results_indexes[key] = sorted(cur_indexes)

    return results, results_indexes

def gene_to_peak_distance(gene_list, all_gtf, peaks_path, bin=1000):
    peaks_dict, peaks_indexes = df_to_index_danpos(peaks_path, bin)

    gene_df = all_gtf[all_gtf['mm9.kgXref.geneSymbol'].isin(gene_list)].copy()
    # print gene_df.shape
    results = []
    for gene in gene_list:
        cur_df = gene_df[gene_df['mm9.kgXref.geneSymbol']==gene]
        best_distance = float('inf')
        for i in range(cur_df.shape[0]):
            if cur_df.iloc[i, 2] == '+':
                cur_tss = cur_df.iloc[i, 3]
            elif cur_df.iloc[i, 2] == '-':
                cur_tss = cur_df.iloc[i, 4]
            cur_chr = cur_df.iloc[i, 1]

            if cur_chr not in peaks_indexes or cur_chr not in peaks_dict:
                continue

            cur_index = bisect_left(peaks_indexes[cur_chr], cur_tss/bin)
            cur_index = cur_index if cur_index < len(peaks_indexes[cur_chr]) else len(peaks_indexes[cur_chr]) -1
            cur_left_index = cur_index - 1 if cur_index - 1 >=0 else 0
            cur_right_index = cur_index + 1 if cur_index + 1 < len(peaks_indexes[cur_chr]) else len(peaks_indexes[cur_chr]) -1

            cur_candidates = peaks_dict[cur_chr][peaks_indexes[cur_chr][cur_index]]
            cur_candidates = cur_candidates.union(peaks_dict[cur_chr][peaks_indexes[cur_chr][cur_left_index]])
            cur_candidates = cur_candidates.union(peaks_dict[cur_chr][peaks_indexes[cur_chr][cur_right_index]])

            # print len(cur_candidates)

            # print gene, best_distance
            for t in cur_candidates:
                cur_dis = abs(cur_tss - t[1])
                # print gene, cur_dis, cur_dis < best_distance
                if cur_dis < best_distance:
                    best_distance = cur_dis
                cur_dis = abs(cur_tss - t[2])
                # print cur_dis
                if cur_dis < best_distance:
                    best_distance = cur_dis
            # print gene, best_distance
        results.append([gene, best_distance])

    df_results = pd.DataFrame(results)
    df_results.columns = ['gene', 'distance']
    return df_results

def saturation_distance(df, name, step=100):
    max_distance = 100000
    total = df.shape[0]
    # print total
    results = []
    for d in range(0, max_distance, step):
        # print d
        cur_df = df[df['distance'] < d]
        results.append([d, cur_df.shape[0]*1.0/total])
    # print results
    results_df = pd.DataFrame(results)
    results_df.to_csv(name, sep='\t', index=None)
    return

all_gtf = pd.read_csv('/archive/tmhkxc48/ref_data/mm9/mm9.20150218.knownGene.xls', sep='\t')
# print all_gtf.columns

gene_list = pd.read_csv('/home/tmhbxx3/archive/JQ1/GSE84520_macrophage/chipseq/wigs/M0_up_list.tsv', sep='\t')
gene_list = set(gene_list.iloc[:, -1].values)

# print gene_list

df1 = gene_to_peak_distance(gene_list, all_gtf, 'SRR049256.Fnor.peaks.xls')

# print df

saturation_distance(df1, 'M0_distance_sat.tsv')

gene_list = pd.read_csv('/home/tmhbxx3/archive/JQ1/GSE84520_macrophage/chipseq/wigs/M1_up_list.tsv', sep='\t')
gene_list = set(gene_list.iloc[:, -1].values)

df2 = gene_to_peak_distance(gene_list, all_gtf, 'SRR049256.Fnor.peaks.xls')

saturation_distance(df2, 'M1_distance_sat.tsv')

gene_list = pd.read_csv('/home/tmhbxx3/archive/JQ1/GSE84520_macrophage/chipseq/wigs/M2_up_list.tsv', sep='\t')
gene_list = set(gene_list.iloc[:, -1].values)

df3 = gene_to_peak_distance(gene_list, all_gtf, 'SRR049256.Fnor.peaks.xls')

saturation_distance(df3, 'M2_distance_sat.tsv')
















