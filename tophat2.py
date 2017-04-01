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
        if len(pair) == 1:
            runTopHatSingle(pair[0], tophatIndex, node_id)
        elif len(pair) == 2:
            runTopHatPair(pair, tophatIndex, node_id)
        finishedjob.append(pair)
        break
    return finishedjob

def runTopHatSingle(sample, tophatIndex, node_id):
    """
    :param sample: sample file path and name
    :param tophatIndex: tophat index path
    :return:
    """
    cmd = "tophat2 -p 8 --mate-std-dev 200 -r 200 "
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

    common_name = find_common(sample1, sample2)

    cmd = "tophat2 -p 8 --mate-std-dev 200 -r 200 "
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
    match = SequenceMatcher(None, string1, string2).find_longest_match(0, len(string1), 0, len(string2))
    if match.a != 0 or match.b != 0:
        return ''
    return string1[match.a:match.a+match.size]

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
        if candidate == string1:
            continue
        cur_score = len(find_common(string1, candidate))
        if cur_score > score:
            match = candidate
            score = cur_score
    return match

def pair_files(path='./'):
    files = [x for x in os.listdir(path) if x.endswith('.fastq') or x.endswith('.gz')]

    pairs = []
    for name in files:
        if name.find('R2') != -1:
            continue
        pairs.append((name, find_most_similar(name, files)))
    return pairs

pairs = pair_files()
print pairs
runTopHat2(pairs)
