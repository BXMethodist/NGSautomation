import os
import pandas as pd


def SRR_single(sample_id, node_id, species='Homo sapien'):
    fastq_dump_cmd = "/home/tmhbxx3/archive/tools/sratoolkit/bin/fastq-dump "
    if species == 'Homo sapiens':
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
    elif species == "Mus musculus":
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/mm9/bowtie/mm9 '

    pbs = open(sample_id+".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N bowtie_"+sample_id+'\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l nodes=1:ppn=8\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd "+os.getcwd()+"\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("wget https://sra-download.ncbi.nlm.nih.gov/srapub/"+sample_id+'\n')
    pbs.write(fastq_dump_cmd + sample_id + '\n')
    pbs.write(bowtie_cmd + sample_id + ".fastq " + sample_id + ".bowtie\n")
    #pbs.write("rm " + bowtie_path + name + ".fastq\n")
    #pbs.write("rm " + bowtie_path + name + "_1.fastq\n")
    #pbs.write("rm " + bowtie_path + name + "_2.fastq\n")
    pbs.close()
    os.system('qsub '+sample_id+".pbs")
    # break
    return

def SRR_pair(sample_id, node_id, species='Homo sapien'):
    fastq_dump_cmd = "/home/tmhbxx3/archive/tools/sratoolkit/bin/fastq-dump --split-3 "
    if species == 'Homo sapiens':
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
    elif species == "Mus musculus":
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/mm9/bowtie/mm9 '

    pbs = open(sample_id + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N bowtie_" + sample_id + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l nodes=1:ppn=8\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("wget https://sra-download.ncbi.nlm.nih.gov/srapub/" + sample_id + '\n')
    pbs.write(fastq_dump_cmd + sample_id + '\n')
    pbs.write('cat ' + sample_id + '_1.fastq ' + sample_id + '_2.fastq >' + sample_id + '.fastq'+'\n')
    pbs.write(bowtie_cmd + sample_id + ".fastq " + sample_id + ".bowtie\n")
    # pbs.write("rm " + bowtie_path + name + ".fastq\n")
    # pbs.write("rm " + bowtie_path + name + "_1.fastq\n")
    # pbs.write("rm " + bowtie_path + name + "_2.fastq\n")
    pbs.close()
    os.system('qsub '+sample_id+".pbs")
    # break
    return

def ENC_pair(sample_id, pair_id, node_id, species='Homo sapien'):
    if species == 'Homo sapiens':
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
    elif species == "Mus musculus":
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/mm9/bowtie/mm9 '
    else:
        print species

    pbs = open(sample_id + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N bowtie_" + sample_id + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l nodes=1:ppn=8\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("wget https://www.encodeproject.org/files/" + sample_id + "/@@download/" + sample_id + ".fastq.gz\n")
    pbs.write("wget https://www.encodeproject.org/files/" + pair_id + "/@@download/" + pair_id + ".fastq.gz\n")
    pbs.write("gunzip "+sample_id + ".fastq.gz\n")
    pbs.write("gunzip " + pair_id + ".fastq.gz\n")
    pbs.write('cat ' + sample_id + '.fastq ' + pair_id + '.fastq >' + sample_id + '.fastq' + '\n')
    pbs.write(bowtie_cmd + sample_id + ".fastq " + sample_id + ".bowtie\n")
    # pbs.write("rm " + bowtie_path + name + ".fastq\n")
    # pbs.write("rm " + bowtie_path + name + "_1.fastq\n")
    # pbs.write("rm " + bowtie_path + name + "_2.fastq\n")
    pbs.close()
    os.system('qsub '+sample_id+".pbs")
    # break
    return

def ENC_single(sample_id, node_id, species='Homo sapien'):
    if species == 'Homo sapiens':
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
    elif species == "Mus musculus":
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/mm9/bowtie/mm9 '

    pbs = open(sample_id + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N bowtie_" + sample_id + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l nodes=1:ppn=8\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("wget https://www.encodeproject.org/files/" + sample_id + "/@@download/" + sample_id + ".fastq.gz\n")
    pbs.write("gunzip "+sample_id + ".fastq.gz\n")
    pbs.write(bowtie_cmd + sample_id + ".fastq " + sample_id + ".bowtie\n")
    # pbs.write("rm " + bowtie_path + name + ".fastq\n")
    # pbs.write("rm " + bowtie_path + name + "_1.fastq\n")
    # pbs.write("rm " + bowtie_path + name + "_2.fastq\n")
    pbs.close()
    os.system('qsub '+sample_id+".pbs")
    # break
    return


def bowtie(search_list, metadata_list):
    """
    run the bowtie for each sample and get the sample input pair list.
    :param search_list:
    :param metadata_list:
    :return: a list of tuples containing pair of sample_id, input_id
    """
    search = pd.read_csv(search_list, index_col=None)
    metadata = pd.read_csv(metadata_list, sep='\t', index_col=None)

    enc_metadata = metadata.set_index(['Run_ID'], drop=False)

    gsm_metadata = metadata.set_index(['GSM_ID'], drop=False)

    search = search.set_index(['Data_ID'], drop=False)

    results = []

    samples = set()

    nodes = [1, 2, 3, 4, 5]

    for i in range(search.shape[0]):
        sample_id = search.ix[i, 'Data_ID']
        input_id = search.ix[i, 'Input']

        if sample_id.startswith('ENC') and not pd.isnull(input_id):
            input_id = input_id[7:-1]

        node_index = i % 5

        species = search.ix[sample_id, 'Organism']

        if sample_id in samples:
            print sample_id
            continue
        else:
            samples.add(sample_id)
            if sample_id.startswith('ENC'):
                if enc_metadata.ix[sample_id, 'Run type'] == 'single-ended':
                    ENC_single(sample_id, nodes[node_index], species)
                elif enc_metadata.ix[sample_id, 'Run type'] == 'paired-ended':
                    pair_id = enc_metadata.ix[sample_id, 'Paired with']
                    ENC_pair(sample_id, pair_id, nodes[node_index], species)
            if sample_id.startswith("GSM"):
                SRR_id = gsm_metadata.ix[sample_id, 'Run_ID']
                if gsm_metadata.ix[sample_id, 'Run type'] == 'SINGLE':
                    SRR_single(SRR_id, node_id=nodes[node_index], species=search.ix[sample_id, 'Organism'])
                elif gsm_metadata.ix[sample_id, 'Run type'] == 'PAIRED':
                    SRR_pair(SRR_id, nodes[node_index], species)


        if pd.isnull(input_id):
            results.append((sample_id, ''))
        else:
            results.append((sample_id, input_id))
            if input_id in samples:
                continue
            else:
                samples.add(input_id)

    df = pd.DataFrame(results)
    df.to_csv('sample_input_pair.csv', index=None, header=False)
    return results

bowtie("Search_ResultHomo_sapiensWithBMI1.csv", "BMI1_chipseq_metadata.txt")


