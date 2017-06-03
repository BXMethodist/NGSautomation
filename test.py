import pandas as pd

# df = pd.read_csv('/Users/boxia/Desktop/reference_map/QC/bowtie_results.csv', index_col=0)
#
# new = []
#
# for i in df.index:
#     if i.find('_') != -1:
#         i = '_'.join(sorted(i.split('_'))).strip()
#         new.append(i)
#     else:
#         new.append(i)
#
#
# print new
# df['new'] = new
# df.set_index(['new'], drop=True)
#
# df.to_csv('sorted_bowtie_results.csv')

# r_df = pd.read_csv('/Users/boxia/Desktop/reference_map/QC/results_NRF.csv', index_col=0)
#
# new = []
#
# for i in r_df.index:
#     if i.find('_') != -1:
#         i = '_'.join(sorted(i.split('_'))).strip()
#         new.append(i)
#     else:
#         new.append(i)
#
# r_df['new2'] = new
# r_df.set_index(['new2'], drop=True)
#
# r_df.to_csv('sorted_results_NRF.csv')
df = pd.read_csv('./sorted_bowtie_results.csv', index_col=0)

r_df = pd.read_csv('./sorted_results_NRF.csv', index_col=0)


result_df = r_df.join(df)

result_df.to_csv("NRF.csv")