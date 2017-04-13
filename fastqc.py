# This module is used to check the quality of fastq file

import os, numpy as np, pandas as pd
from collections import defaultdict

def RandomReads(fastq, number_reads, individual):
    """
    Generate number of random reads to to the NCBI blastn
    :param fastq: the fastq file path
    :param number_reads: random number of reads
    :return:
    """

    fastq_obj = open(fastq, 'r')
    info = fastq_obj.readlines()
    fastq_obj.close()

    total_reads = len(info)/4

    candidates = list(np.random.randint(10000, total_reads-10000, number_reads))

    results = []

    for candidate in candidates:
        start = candidate * 4
        end = candidate * 4 + 4
        reads = info[start: end]

        results.append('>read'+str(candidate))
        results.append(reads[1])
    name = fastq[:-6]
    if individual:
        results_obj = open(name+'_random_reads.txt', 'w')
        for r in results:
            results_obj.write(r+'\n')
        results_obj.close()
    return results

def BatchRandomReads(fastqs, number_reads, surffix, individual=False):
    """
    Generate a batch of fastq files' random number of reads
    :param fastqs: a folder containing fastq files
    :param number_reads: number of random reads
    :param file name identifier
    :return:
    """
    files = [x for x in os.listdir(fastqs) if x.endswith('.fastq') and x.find(surffix)!= -1]
    results = []
    for f in files:
        result = RandomReads(f, number_reads, individual)
        results += result

    results_obj = open('total' + '_random_reads.txt', 'w')
    for r in results:
        results_obj.write(r + '\n')
    results_obj.close()
    return results

def RunFastqc(fastqs):
    os.system('module load fastqc')
    files = [x for x in os.listdir(fastqs) if x.endswith('.fastq')]
    os.system('fastqc '+ ' '.join(files))

def Blast_Organism(file, count=1):
    """
    :param file: output txt file from blastn
    :param count: top n species from blastn results
    :return:
    """
    f = open(file, 'r')
    info = f.readlines()
    f.close()
    species = defaultdict(int)
    genuses = defaultdict(int)
    cur_count = 0
    for line in info:
        if line.startswith('ALIGNMENTS'):
            cur_count = 0
        if line.startswith('>') and cur_count < count:
            line = line[1:].split()
            if line[1] == "PREDICTED:":
                species_name = line[2] + " " + line[3]
                genus_name = line[2]
            else:
                species_name = line[1] + " "+ line[2]
                genus_name = line[1]
            species[species_name] += 1
            genuses[genus_name] +=1
            cur_count += 1

    df = pd.DataFrame.from_dict(species, orient='index')
    df.columns = ['count']
    df.sort(columns=['count'], inplace=True, ascending=False)

    df.to_csv('pancrease_species_for_reads.csv')

    df = pd.DataFrame.from_dict(genuses, orient='index')
    df.columns = ['count']
    df.sort(columns=['count'], inplace=True, ascending=False)

    df.to_csv('pancrease_genuses_for_reads.csv')
    return df

Blast_Organism('/Users/boxia/Desktop/TERT/pancrease_reads_blast.txt')





# BatchRandomReads("./", 100, 'P')

#
# df = pd.read_csv('~/Desktop/TERT/tophat2_result.csv', index_col=0)
#
# print np.corrcoef(df, rowvar = False)
#
# print df.columns
# df = Blast_Organism('ES19PKYT016-Alignment.txt')
#
# df = df.ix[:5, :]
#
# import matplotlib.pyplot as plt
#
# df.plot.pie(y='count')
# plt.legend().draggable()
# plt.axis('off')
# plt.show()







