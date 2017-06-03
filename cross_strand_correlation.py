import os, pandas as pd

df = pd.read_csv('ENC_H3K4me3_sample_pairs.csv', index_col=0)

df = df[pd.notnull(df['input'])]

tagAligns = [index + '.tagAlign' for index in df.index]
path = '/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bed/'

print len(tagAligns)
# print bowties[107]

for tagAlign in tagAligns[:200]:
    os.system('cp '+path+tagAlign + ' ./')
    cmd = 'Rscript /share/apps/phantompeakqualtools/run_spp.R -c=/home/tmhbxx3/archive/H3K4me3/ENCODE_with_input/NSC_RSC/' \
          + tagAlign +' -savp -out=' + tagAlign[:-9]+'.txt'

    pbs = open(tagAlign + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N NSC/RSC_" + tagAlign[:-9] + '\n')
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l pmem=4000mb\n")
    # pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("module load R/3.1.2\n")
    pbs.write("module load phantompeakqualtools\n")
    pbs.write(cmd + '\n')
    pbs.close()
    os.system('qsub ' + tagAlign + ".pbs")

    # break



# Rscript /share/apps/phantompeakqualtools/run_spp.R -c=/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bed/ENCFF387HNJ.tagAlign -savp -out=result.txt