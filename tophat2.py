import os
import pandas as pd
import numpy as np
from difflib import SequenceMatcher

def runTopHat2(samplelists, finishedjob=[], tophatIndex=" /archive/tmhkxc48/ref_data/mm9/bowtie2/mm9 "):
    """
    :param samplelists: a list of tuples, if tuple length is one then single, 2 then pair
    :param finishedjob: list of finished jobs
    :param tophatIndex: tophat index path
    :return:
    """
    nodes = [1,2,3,4,5, 6]
    n = 0
    for pair in samplelists:
        node_id = nodes[n%6]
        if pair in finishedjob:
            continue
        n += 1
        if len(pair) == 1 or pair[1] is None:
            runTopHatSingle(pair[0], tophatIndex, node_id)
        elif len(pair) == 2:
            runTopHatPair(pair, tophatIndex, node_id)
        finishedjob.append(pair)
        #break
    return finishedjob

def runTopHatSingle(sample, tophatIndex, node_id):
    """
    :param sample: sample file path and name
    :param tophatIndex: tophat index path
    :return:
    """
    cmd = "tophat2 -p 8 "
    if sample.find("/") == -1:
        outputstartindex = 0
    else:
        outputstartindex = sample.rfind("/") + 1
    outputname = sample[outputstartindex:sample[outputstartindex:].find('.')]

    if not os.path.isdir('./'+outputname):
        os.system("mkdir "+'./'+outputname)

    cmd += '-o '+'./'+outputname +' ' + tophatIndex + ' ' + sample

    submit_pbs(cmd, outputname, node_id)


def runTopHatPair(pair, tophatIndex, node_id):
    """
    :param pair: pair end sample, a tuple have two elements
    :param tophatIndex: tophat index path
    :return:
    """
    sample1, sample2 = pair

    common_name = find_common_name(sample1, sample2)

    cmd = "tophat2 -p 8 "
    if common_name.find("/") == -1:
        outputstartindex = 0
    else:
        outputstartindex = common_name.rfind("/") + 1
    outputname = common_name[outputstartindex:common_name[outputstartindex:].find('.')]

    if not os.path.isdir('./' + outputname):
        os.system("mkdir " + './' + outputname)

    cmd += '-o ' + './' + outputname + ' ' + tophatIndex + ' ' + sample1 +' ' +  sample2

    submit_pbs(cmd, outputname, node_id)

def submit_pbs(cmd, outputname, node_id):
    """
    :param cmd: command line
    :param outputname: output folder name
    :return:
    """
    pbs = open(outputname + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N " + outputname + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l nodes=1:ppn=8\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write(cmd + "\n")
    pbs.close()
    os.system('qsub ' + outputname + ".pbs")

def find_common(string1, string2):
    """
    :param name1: string 1
    :param name2: string 2
    :return: longest common part
    """
    if string2.find('_1.fq.gz') != -1:
        return 0

    s1 = string1.replace('_1.fq.gz', '')
    s2 = string2.replace('_2.fq.gz', '')
    # print s1, s2
    match = SequenceMatcher(None, s1, s2).find_longest_match(0, len(s1), 0, len(s2))
    # print match.size
    return match.size
    # print match.size, match.a, match.b
    # if match.a != 0 or match.b != 0:
    #     return 0
    # return match.size

def find_common_name(string1, string2):
    """
    :param name1: string 1
    :param name2: string 2
    :return: longest common part
    """
    if string2.find('R1') != -1:
        return 0

    s1 = string1.replace('_1.fq.gz', '')
    s2 = string2.replace('_2.fq.gz', '')
    match = SequenceMatcher(None, s1, s2).find_longest_match(0, len(s1), 0, len(s2))
    # if match.a != 0 or match.b != 0:
    #     return 0
    return s1[match.a:match.a+match.size]

def find_most_similar(string1, strings):
    """
    find the most similar string with the longest common part
    :param string1:
    :param strings:
    :return: a most similar string
    """
    score = 0
    match = None
    for candidate in strings:
        if candidate == string1 or candidate.find('_1.fq.gz')!= -1:
            continue
        # print string1, candidate
        cur_score = find_common(string1, candidate)
        if cur_score > score:
            match = candidate
            score = cur_score
    return match

def pair_files(path='./'):
    files = [x for x in os.listdir(path) if x.endswith('.fastq') or x.endswith('.gz')]

    pairs = []
    for name in files:
        if name.find('_2.fq.gz') != -1:
            continue
        pairs.append((name, find_most_similar(name, files)))
    return pairs

def tophat2_single_result(directories):
    """
    :param directories: list of directories
    :return:
    """
    columns = ['total_reads', 'mapped_reads']
    results = []
    folders_names = []
    for folder in directories:
        folder = folder[:-1] if folder.endswith('/') else folder
        try:
            file_name = folder + "/align_summary.txt"
            f = open(file_name, "r")
            info = f.readlines()
            f.close()
        except:
            continue

        name_index = folder.rfind('/') + 1 if folder.rfind('/') != -1 else 0
        folder_name = folder[name_index:]
        folders_names.append(folder_name)

        cur_total_reads = info[1][info[1].find(":") + 1:].strip()
        cur_mapped_reads = info[2][info[2].find(":") + 1:info[2].find("(")].strip()
        results.append((cur_total_reads, cur_mapped_reads))
    df = pd.DataFrame(results, index=folders_names, columns=columns)

    df.to_csv("tophat2_single_result.csv")
    return df

def tophat2_pair_result(directories):
    """
    :param directories: list of directories
    :return:
    """
    columns = ['total_reads', 'mapped_reads']
    results = []
    folders_names = []
    for folder in directories:
        folder = folder[:-1] if folder.endswith('/') else folder
        try:
            file_name = folder + "/align_summary.txt"
            f = open(file_name, "r")
            info = f.readlines()
            f.close()
        except:
            continue

        name_index = folder.rfind('/') + 1 if folder.rfind('/') != -1 else 0
        folder_name = folder[name_index:]
        folders_names.append(folder_name)

        cur_total_reads = info[1][info[1].find(":")+1:].strip()
        cur_mapped_reads = info[10][info[10].find(":")+1:].strip()
        results.append((cur_total_reads, cur_mapped_reads))
    df = pd.DataFrame(results, index=folders_names, columns=columns)

    df.to_csv("tophat2_pair_result.csv")
    return df

def moveAndChangeNameForTophat(path, target_path="../bams"):
    '''
    move bam files from folders and transfer bam to the target_path
    :param path: the directory contains tophat_output
    :param target_path:
    :return:
    '''
    if not os.path.isdir(target_path):
        os.system("mkdir " + target_path)

    for folder in path:
        old_file_location = folder+"/accepted_hits.bam"
        new_file_location = target_path+"/"+folder+".bam"
        os.system("mv " + old_file_location + " "+ new_file_location)
    return

# directories = [x for x in os.listdir("./") if os.path.isdir(x) and x.endswith("001")]
#
# moveAndChangeNameForTophat([x for x in os.listdir('.') if x.endswith('_')])
#
# pairs = pair_files('../FASTQ')
# for p in pairs:
#     print p

# print pair_files()

runTopHat2(pair_files())

# tophat2_pair_result(os.listdir('.'))