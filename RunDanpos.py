# This module contains functions to run danpos and related function
import os, pandas as pd

def danpos_no_input(sample_id, search_df, metadata_df, node_id):
    if sample_id.startswith("GSM"):
        sample_id = metadata_df.ix[sample_id, 'Run_ID']

    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dregion '
    danpos_parameters = ' -u 1 --smooth_width 0 -c 25000000 --frsz 200 --extend 200 ' \
                        '--extend_dis 3000 --pheight 1e-8 -ep 1e-5 -o ' + os.getcwd() + '/' + sample_id

    cmd = danpos_cmd + sample_id+'.bowtie' + danpos_parameters

    pbs = open(sample_id + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N danpos_" + sample_id + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    # pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("module load R/3.2.1\n")
    pbs.write(cmd + '\n')
    pbs.close()
    os.system('qsub ' + sample_id + ".pbs")
    return

def danpos_input(sample_id, input_id, search_df, metadata_df, node_id):
    if sample_id.startswith("GSM"):
        sample_id = metadata_df.ix[sample_id, 'Run_ID']
        input_id = metadata_df.ix[input_id, 'Run_ID']

    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dregion '
    danpos_parameters = ' -u 1 --smooth_width 0 -c 25000000 --frsz 200 --extend 200 ' \
                        '--extend_dis 3000 --pheight 1e-8 -ep 1e-5 -o ' + os.getcwd() + '/' + sample_id

    cmd = danpos_cmd + sample_id + '.bowtie' +' -b '+ input_id +".bowtie" + danpos_parameters

    pbs = open(sample_id + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N danpos_" + sample_id + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    # pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("module load R/3.2.1\n")
    pbs.write(cmd + '\n')
    pbs.close()
    os.system('qsub ' + sample_id + ".pbs")
    return

def RunDanpos(search, metadata, sample_input="sample_input_pair.csv"):
    """
    :param sample_input: file path for sample input pair
    :param search: search dataframe from chipseqpair search function
    :param metadata: metadata dataframe from chipseqpair query function
    :return:
    """
    sample_input_df = pd.read_csv(sample_input, header=None, index_col=None)
    search_df = pd.read_csv(search, index_col=None)
    metadata_df = pd.read_csv(metadata, index_col=None, sep='\t')
    metadata_df = metadata_df.set_index(['GSM_ID'])

    nodes = [1,2,3,4,5, 6]

    for i in range(sample_input_df.shape[0]):
        node_id = nodes[i%6]
        sample_id = sample_input_df.ix[i, 0]
        input_id = sample_input_df.ix[i, 1]
        if pd.isnull(input_id):
            danpos_no_input(sample_id, search_df, metadata_df, node_id)
        else:
            danpos_input(sample_id, input_id, search_df, metadata_df, node_id)


RunDanpos("Search_ResultHomo_sapiensWithBMI1_ENC.csv", "BMI1_chipseq_metadata.txt")