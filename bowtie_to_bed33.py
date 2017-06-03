import os, pandas as pd

df = pd.read_csv('ENC_H3K4me3_sample_pairs.csv', index_col=0)

df = df[pd.notnull(df['input'])]

bowties = [index + '.bowtie' for index in df.index]
path = '/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/'

print len(bowties)
# print bowties[107]

for bowtie in bowties[:200]:
    f = open(path + bowtie, 'r')
    result = open(bowtie[:-7] +'.tagAlign', 'w')

    for line in f:
        line = [x.strip() for x in line.split('\t')]
        index = None

        if line[1] == '-' or line[1] == '+':
            index = 1
        elif line[2] == '-' or line[1] == '+':
            index = 2
        if index is None:
            print 'something wrong!!!'
        else:
            result.write(line[index+1] + '\t' + line[index+2] + '\t' + str(len(line[index+3])+int(line[index+2])) + '\t' +'N' + '\t' +'111'+'\t'+line[index] +'\n')

    f.close()
    result.close()

    # cmd = 'Rscript /share/apps/phantompeakqualtools/run_spp.R -c=/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bed/' \
    #       + bowtie[:-7]+'.tagAlign' +' -savp -out=result.txt'
    #
    # os.system(cmd)


# Rscript /share/apps/phantompeakqualtools/run_spp.R -c=/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bed/ENCFF387HNJ.tagAlign -savp -out=result.txt