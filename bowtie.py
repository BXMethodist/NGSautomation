import os
import pandas as pd, urllib2, json


def ENC_biosample_meta(biosample_id):
    # print biosample_id
    if biosample_id is None or pd.isnull(biosample_id):
        return None

    url = "https://www.encodeproject.org/experiments/" + biosample_id + "/?format=json"
    try:
        response = urllib2.urlopen(url)
        data = json.loads(response.read())
    except:
        print url
        return None
    return data

def SRR_single(sample_id, node_id, species='Homo sapiens'):
    """
        submit a pbs for SRR single end to generate bowtie result
        :param sample_id: SRR ID
        :param node_id: specify the node, sometimes some node does not work.
        :param species: the species, only for human and mouse
        :return:
        """
    fastq_dump_cmd = "/home/tmhbxx3/archive/tools/sratoolkit/bin/fastq-dump "
    if species == 'Homo sapiens':
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
    elif species == "Mus musculus":
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/mm9/bowtie/mm9 '

    pbs = open(sample_id+".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N bowtie_"+sample_id+'\n')
    pbs.write("#PBS -q highmem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=32:00:00\n")
    pbs.write("#PBS -l nodes=1:ppn=8\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    #pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd "+os.getcwd()+"\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write(fastq_dump_cmd + sample_id + '\n')
    pbs.write(bowtie_cmd + sample_id + ".fastq " + sample_id + ".bowtie\n")
    pbs.write("rm " + sample_id + ".fastq\n")
    pbs.close()
    os.system('qsub '+sample_id+".pbs")
    # break
    return

def SRR_pair(sample_id, node_id, species='Homo sapiens'):
    """
        submit a pbs for SRR pair end to generate bowtie result
        :param sample_id: SRR ID
        :param node_id: specify the node, sometimes some node does not work.
        :param species: the species, only for human and mouse
        :return:
        """
    fastq_dump_cmd = "/home/tmhbxx3/archive/tools/sratoolkit/bin/fastq-dump --split-3 "
    if species == 'Homo sapiens':
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
    elif species == "Mus musculus":
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/mm9/bowtie/mm9 '

    pbs = open(sample_id + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N bowtie_" + sample_id + '\n')
    pbs.write("#PBS -q highmem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=32:00:00\n")
    pbs.write("#PBS -l nodes=1:ppn=8\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    #pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    # pbs.write("wget https://sra-download.ncbi.nlm.nih.gov/srapub/" + sample_id + '\n')
    pbs.write(fastq_dump_cmd + sample_id + '\n')
    pbs.write('cat ' + sample_id + '_1.fastq ' + sample_id + '_2.fastq >' + sample_id + '.fastq'+'\n')
    pbs.write(bowtie_cmd + sample_id + ".fastq " + sample_id + ".bowtie\n")
    pbs.write("rm " + sample_id + ".fastq\n")
    pbs.write("rm " + sample_id + "_1.fastq\n")
    pbs.write("rm " + sample_id + "_2.fastq\n")
    pbs.close()
    os.system('qsub '+sample_id+".pbs")
    # break
    return

def ENC_pair(sample_id, pair_id, node_id, species='Homo sapiens'):
    """
    :param sample_id: ENC ID
    :param pair_id: the corresponding pair_id for pair end sequencing
    :param node_id: specify the node, sometimes some node does not work.
    :param species: the species, only for human and mouse
    :return:
    """
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
    #pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")

    if not os.path.isfile(sample_id+'.fastq'):
        pbs.write("wget https://www.encodeproject.org/files/" + sample_id + "/@@download/" + sample_id + ".fastq.gz\n")
        pbs.write("gunzip " + sample_id + ".fastq.gz\n")
    if not os.path.isfile(pair_id+'.fastq'):
        pbs.write("wget https://www.encodeproject.org/files/" + pair_id + "/@@download/" + pair_id + ".fastq.gz\n")
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

def ENC_single(sample_id, node_id, species='Homo sapiens'):
    """
    submit a pbs for ENCODE single end to generate bowtie result
    :param sample_id: ENC ID
    :param node_id: specify the node, sometimes some node does not work.
    :param species: the species, only for human and mouse
    :return:
    """
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
    #pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")

    if not os.path.isfile(sample_id+'.fastq'):
        pbs.write("wget https://www.encodeproject.org/files/" + sample_id + "/@@download/" + sample_id + ".fastq.gz\n")
        pbs.write("gunzip " + sample_id + ".fastq.gz\n")

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
    :param search_list: search results from ChipSeqpair search function
    :param metadata_list: query results from ChipSeqpair query function
    :return: a list of tuples containing pair of sample_id, input_id
    """
    search = pd.read_csv(search_list, index_col=None)
    metadata = pd.read_csv(metadata_list, sep='\t', index_col=None)

    enc_metadata = metadata.set_index(['Run_ID'], drop=False)

    gsm_metadata = metadata.set_index(['GSM_ID'], drop=False)
    gsm_metadata = gsm_metadata.drop_duplicates(subset=['Run_ID'])

    search = search.set_index(['Data_ID'], drop=False)

    results = []

    samples = set()

    nodes = [1, 2, 3, 4, 5, 6]

    missed =set()
    count = 0

    for i in range(search.shape[0]):
        sample_id = search.ix[i, 'Data_ID']
        input_id = search.ix[i, 'Input']

        # if pd.isnull(input_id):
        #     continue

        if sample_id.startswith('ENC') and not pd.isnull(input_id):
            input_id = input_id[7:-1]

        node_index = i % 6

        species = search.ix[sample_id, 'Organism']

        if sample_id.startswith('ENC'):
            try:
                if enc_metadata.ix[sample_id, 'Run type'] == 'single-ended':
                    if sample_id not in samples:
                        # ENC_single(sample_id, nodes[node_index], species)
                        samples.add(sample_id)
                    else:
                        continue
                elif enc_metadata.ix[sample_id, 'Run type'] == 'paired-ended':
                    pair_id = enc_metadata.ix[sample_id, 'Paired with']
                    if sample_id not in samples and pair_id not in samples:
                        samples.add(sample_id)
                        samples.add(pair_id)
                        # ENC_pair(sample_id, pair_id, nodes[node_index], species)
                    elif sample_id in samples:
                        samples.add(pair_id)
                        continue
                    elif pair_id in samples:
                        samples.add(sample_id)
                        continue
            except:
                print sample_id, 'sample'
                missed.add(sample_id)
        if sample_id.startswith("GSM"):
            if sample_id in samples:
                continue
            samples.add(sample_id)

            try:
                SRR_ids = gsm_metadata.ix[sample_id, 'Run_ID']
                if isinstance(SRR_ids, str):
                    if gsm_metadata.ix[sample_id, 'Run type'] == 'SINGLE':
                        # SRR_single(SRR_id, node_id=nodes[node_index], species=species)
                        pass
                    elif gsm_metadata.ix[sample_id, 'Run type'] == 'PAIRED':
                        # SRR_pair(SRR_id, nodes[node_index], species)
                        pass
                else:
                    for srr_id in SRR_ids.tolist():
                        if gsm_metadata[gsm_metadata.Run_ID == srr_id].ix[0, 'Run type'] == 'SINGLE':
                            # SRR_single(SRR_id, node_id=nodes[node_index], species=species)
                            pass
                        elif gsm_metadata[gsm_metadata.Run_ID == srr_id].ix[0, 'Run type'] == 'PAIRED':
                            # SRR_pair(SRR_id, nodes[node_index], species)
                            pass
            except:
                print sample_id, 'sample'
                missed.add(sample_id)
                pass

        if pd.isnull(input_id):
            # results.append((sample_id, ''))
            pass
        else:
            input_id = input_id.replace('files','')
            input_id = input_id.replace('/','')
            input_ids = [x.strip() for x in input_id.split(',')]
            if len(input_ids) > 1:
                print sample_id
                if sample_id.startswith('ENC'):
                    count +=1
            if sample_id.startswith('ENC'):
                results.append((sample_id, ';'.join(input_ids)))
            if os.path.isfile(input_id):
                continue
            else:
                for input_id in input_ids:
                    samples.add(input_id)
                    if input_id.startswith('ENC'):
                        try:
                            if enc_metadata.ix[input_id, 'Run type'] == 'single-ended':
                                # ENC_single(input_id, nodes[node_index], species)
                                pass
                            elif enc_metadata.ix[input_id, 'Run type'] == 'paired-ended':
                                pair_id = enc_metadata.ix[input_id, 'Paired with']
                                samples.add(pair_id)
                                # ENC_pair(input_id, pair_id, nodes[node_index], species)
                                pass
                        except:
                            print input_id, 'input'
                            missed.add(input_id)
                    if input_id.startswith("GSM"):
                        try:
                            SRR_ids = gsm_metadata.ix[input_id, 'Run_ID']
                            if isinstance(SRR_ids, str):
                                if gsm_metadata.ix[input_id, 'Run type'] == 'SINGLE':
                                    # SRR_single(SRR_id, node_id=nodes[node_index], species=species)
                                    pass
                                elif gsm_metadata.ix[input_id, 'Run type'] == 'PAIRED':
                                    # SRR_pair(SRR_id, nodes[node_index], species)
                                    pass
                            else:
                                for srr_id in SRR_ids.tolist():
                                    if gsm_metadata[gsm_metadata.Run_ID == srr_id].ix[0, 'Run type'] == 'SINGLE':
                                        # SRR_single(SRR_id, node_id=nodes[node_index], species=species)
                                        pass
                                    elif gsm_metadata[gsm_metadata.Run_ID == srr_id].ix[0, 'Run type'] == 'PAIRED':
                                        # SRR_pair(SRR_id, nodes[node_index], species)
                                        pass
                        except:
                            print input_id, 'input'
                            missed.add(input_id)
                            pass


    df = pd.DataFrame(results)
    df.to_csv('sample_input_pair.csv', index=None, header=False)

    f = open('missed.csv', 'w')
    for id in set(missed):
        f.write(id+'\n')
    f.close()
    print count
    return results

def GEO_sample_input_pair(search_list, metadata_list):
    """
    run the bowtie for each sample and get the sample input pair list.
    :param search_list: search results from ChipSeqpair search function
    :param metadata_list: query results from ChipSeqpair query function
    :return: a list of tuples containing pair of sample_id, input_id
    """
    search = pd.read_csv(search_list, index_col=None)
    metadata = pd.read_csv(metadata_list, sep='\t', index_col=None)

    enc_metadata = metadata.set_index(['Run_ID'], drop=False)

    gsm_metadata = metadata.set_index(['GSM_ID'], drop=False)
    gsm_metadata = gsm_metadata.drop_duplicates(subset=['Run_ID'])

    search = search.set_index(['Data_ID'], drop=False)

    results = []

    samples = set()

    nodes = [1, 2, 3, 4, 5, 6]

    missed =set()
    count = 0

    for i in range(search.shape[0]):
        sample_id = search.ix[i, 'Data_ID']
        input_id = search.ix[i, 'Input']

        if pd.isnull(input_id):
            continue

        if sample_id.startswith('ENC'):
            continue

        node_index = i % 6

        species = search.ix[sample_id, 'Organism']

        if sample_id.startswith("GSM"):
            if sample_id in samples:
                continue
            samples.add(sample_id)

            try:
                SRR_ids = gsm_metadata.ix[sample_id, 'Run_ID']
                if isinstance(SRR_ids, str):
                    if gsm_metadata.ix[sample_id, 'Run type'] == 'SINGLE':
                        # SRR_single(SRR_id, node_id=nodes[node_index], species=species)
                        pass
                    elif gsm_metadata.ix[sample_id, 'Run type'] == 'PAIRED':
                        # SRR_pair(SRR_id, nodes[node_index], species)
                        pass
                else:
                    for srr_id in SRR_ids.tolist():
                        if gsm_metadata[gsm_metadata.Run_ID == srr_id].ix[0, 'Run type'] == 'SINGLE':
                            # SRR_single(SRR_id, node_id=nodes[node_index], species=species)
                            pass
                        elif gsm_metadata[gsm_metadata.Run_ID == srr_id].ix[0, 'Run type'] == 'PAIRED':
                            # SRR_pair(SRR_id, nodes[node_index], species)
                            pass
            except:
                print sample_id, 'sample'
                missed.add(sample_id)
                pass

        input_ids = [x.strip() for x in input_id.split(',')]
        if len(input_ids) > 1:
            print sample_id, input_ids
            results.append((sample_id, ';'.join(input_ids)))
        else:
            results.append((sample_id, input_ids[0]))

        if os.path.isfile(input_id):
            continue
        else:
            for input_id in input_ids:
                samples.add(input_id)
                if input_id.startswith('ENC'):
                    try:
                        if enc_metadata.ix[input_id, 'Run type'] == 'single-ended':
                            # ENC_single(input_id, nodes[node_index], species)
                            pass
                        elif enc_metadata.ix[input_id, 'Run type'] == 'paired-ended':
                            pair_id = enc_metadata.ix[input_id, 'Paired with']
                            samples.add(pair_id)
                            # ENC_pair(input_id, pair_id, nodes[node_index], species)
                            pass
                    except:
                        print input_id, 'input'
                        missed.add(input_id)
                if input_id.startswith("GSM"):
                    try:
                        SRR_ids = gsm_metadata.ix[input_id, 'Run_ID']
                        if isinstance(SRR_ids, str):
                            if gsm_metadata.ix[input_id, 'Run type'] == 'SINGLE':
                                # SRR_single(SRR_id, node_id=nodes[node_index], species=species)
                                pass
                            elif gsm_metadata.ix[input_id, 'Run type'] == 'PAIRED':
                                # SRR_pair(SRR_id, nodes[node_index], species)
                                pass
                        else:
                            for srr_id in SRR_ids.tolist():
                                if gsm_metadata[gsm_metadata.Run_ID == srr_id].ix[0, 'Run type'] == 'SINGLE':
                                    # SRR_single(SRR_id, node_id=nodes[node_index], species=species)
                                    pass
                                elif gsm_metadata[gsm_metadata.Run_ID == srr_id].ix[0, 'Run type'] == 'PAIRED':
                                    # SRR_pair(SRR_id, nodes[node_index], species)
                                    pass
                    except:
                        print input_id, 'input'
                        missed.add(input_id)
                        pass


    df = pd.DataFrame(results)
    df.to_csv('sample_input_pair.csv', index=None, header=False)

    f = open('missed.csv', 'w')
    for id in set(missed):
        f.write(id+'\n')
    f.close()
    print count
    return results

def ENC_multiple_bowtie(ids, species='Homo sapiens'):
    if species == 'Homo sapiens':
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
    elif species == "Mus musculus":
        bowtie_cmd = 'bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/mm9/bowtie/mm9 '
    else:
        print species

    run_name = '_'.join(ids)
    pbs = open(run_name + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N " + run_name + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l nodes=1:ppn=8\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    #pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")

    for id in ids:
        if not os.path.isfile(id+'.fastq'):
            pbs.write("wget -q https://www.encodeproject.org/files/" + id + "/@@download/" + id + ".fastq.gz\n")
            pbs.write("gunzip " + id + ".fastq.gz\n")

    if len(ids) > 1:
        cat_cmd = ' '.join([id + '.fastq' for id in ids])
        if not os.path.isfile(run_name + '.fastq'):
            pbs.write('cat ' + cat_cmd + ' >' + run_name + '.fastq' + '\n')
    pbs.write(bowtie_cmd + run_name + ".fastq " + run_name + ".bowtie\n")
    pbs.close()
    os.system('qsub '+run_name+".pbs")
    return


def bowtie_ENC(ENC_metadata_list):
    """
    :param ENC_metadata_list: metadata from ENC
    :return:
    """
    samples = {}
    bam_not_covered_exp = set()
    no_derived_from = set()
    no_files = set()


    # step 1 group the target dataset by experiment ID and biological replicates ID, biological pairs

    # url = 'https://www.encodeproject.org/metadata/type=Experiment&files.file_type=fastq/metadata.tsv'
    # df = pd.read_csv(url, sep='\t', dtype=str, index_col=0)

    df = pd.read_csv('all_ENC_metadata.csv', index_col=0)

    # bam_url = 'https://www.encodeproject.org/metadata/type=Experiment&files.file_type=bam/metadata.tsv'
    # bam_df = pd.read_csv(bam_url, sep='\t', dtype=str, index_col=0)

    sample_df = df[(df['Experiment target'] == 'H3K4me3-human') & (df['Assay'] == 'ChIP-seq')]

    # print sample_df.index

    #     exp_metadata = ENC_biosample_meta(experiment)
    #     bam_count = 0
    #     if 'files' not in exp_metadata:
    #         no_files.add((experiment, 'sample'))
    #     for alignment in exp_metadata['files']:
    #         if alignment['file_type'] != 'bam' and alignment['output_type'] != 'unfiltered alignments' and \
    #                 alignment['output_category'] != 'alignment':
    #             continue
    #
    #         if 'derived_from' not in alignment:
    #             print "no derived from information?", alignment['accession'], experiment
    #             no_derived_from.add((alignment['accession'], experiment, 'sample'))
    #             continue
    #
    #         cur_fastqs = [x.replace('/files/', '').replace('/','').strip() for x in alignment['derived_from']]
    #         print cur_fastqs
    #
    #
    #         final_input_fastqs = set()
    #
    #         related_input_exp = set()
    #
    #         related_input_fastq = set()
    #
    #         for fastq in cur_fastqs:
    #             if fastq not in sample_df.index:
    #                 continue
    #             if not pd.isnull(sample_df.ix[fastq, 'Controlled by']):
    #                 known_inputs = [input.replace('/', '').replace('files', '') for input in
    #                                 sample_df.ix[fastq, 'Controlled by'].split(',')]
    #
    #                 for cur_input in known_inputs:
    #                     related_input_fastq.add(cur_input)
    #                     if cur_input in df.index:
    #                         if not pd.isnull(df.ix[cur_input, 'Paired with']):
    #                             paired_input = df.ix[cur_input, 'Paired with']
    #                             related_input_fastq.add(paired_input.strip())
    #                             if paired_input in df.index:
    #                                 related_input_exp.add(df.ix[paired_input, 'Experiment accession'])
    #                         related_input_exp.add(df.ix[cur_input, 'Experiment accession'])
    #         for input_exp in related_input_exp:
    #             cur_input_meta = ENC_biosample_meta(input_exp)
    #             if 'files' not in cur_input_meta:
    #                 print 'no files in metada?!', input_exp
    #                 no_files.add((input_exp, 'input'))
    #                 continue
    #             for input_alignment in cur_input_meta['files']:
    #                 if input_alignment['file_type'] != 'bam' and input_alignment['output_type'] != 'unfiltered alignments' and \
    #                                 input_alignment['output_category'] != 'alignment':
    #                     continue
    #
    #                 if 'derived_from' not in input_alignment:
    #                     print "no derived from information?", input_alignment['accession'], input_exp
    #                     no_derived_from.add((input_alignment['accession'], input_exp, 'input'))
    #                     continue
    #                 cur_input_fastqs = [x.replace('/files/', '').replace('/','').strip() for x in input_alignment['derived_from']]
    #
    #                 input_real = False
    #                 for cur_input_fastq in cur_input_fastqs:
    #                     if cur_input_fastq in related_input_fastq:
    #                         input_real = True
    #                         break
    #                 if input_real:
    #                     final_input_fastqs.add(tuple(cur_input_fastqs))
    #         samples[tuple(cur_fastqs)] = final_input_fastqs
    #
    #         # cur_samples.add(sample_df.ix[id, 'Paired with'])
    #
    #         bam_count += 1
    #     if bam_count == 0:
    #         print "no bam ?!", experiment
    #         bam_not_covered_exp.add(experiment)
    #
    #
    # print len(samples)
    # print 'no bam', bam_not_covered_exp
    # print 'no files', no_files
    # print 'no derived from', no_derived_from

    for experiment in sample_df['Experiment accession'].unique():
        replicates = sample_df[sample_df['Experiment accession'] == experiment]['Biological replicate(s)'].unique()
        for rep in replicates:
            # try:
            cur_samples = sample_df[(sample_df['Experiment accession'] == experiment) & (sample_df['Biological replicate(s)'] == rep)].index
            cur_samples = set(cur_samples)

            known_samples = [x.strip() for x in cur_samples]

            for id in known_samples:
                if not pd.isnull(sample_df.ix[id, 'Paired with']):
                    cur_samples.add(sample_df.ix[id, 'Paired with'].strip())

            # if len(sample_df.ix[list(cur_samples), 'Library ID'].unique()) != 1:
            #     print cur_samples
            #     print sample_df.ix[list(cur_samples), ['Library ID', 'Biological Replicate ID', 'Technical Replicate ID']]

            cur_inputs = set()

            for id in cur_samples:
                if not pd.isnull(df.ix[id, 'Controlled by']):
                    known_inputs = [input.replace('/','').replace('files', '') for input in df.ix[id, 'Controlled by'].split(',')]

                    for known_id in known_inputs:
                        known_id = known_id.strip()
                        # if known_id == "ENCFF377FRL":
                        #     print pd.isnull(df.ix[known_id, 'Paired with']), df.ix[known_id, 'Paired with']
                        # print known_id
                        cur_inputs.add(known_id.strip())
                        if known_id not in df.index:
                            continue
                        if not pd.isnull(df.ix[known_id, 'Paired with']):
                            cur_inputs.add(df.ix[known_id, 'Paired with'].strip())
                else:
                    pass

            cur_samples = [x.strip() for x in cur_samples]
            # cur_inputs = [tuple(y) for y in cur_inputs]



            samples[tuple(set(cur_samples))] = tuple(cur_inputs)
            # except:
            #     continue
    # step 2 cat the fastq in each group and run bowtie
    results = []

    count = 0
    finished = [x for x in os.listdir('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/') if x.endswith('bowtie')]
    # finished = []
    for i in range(len(finished)):
        x = finished[i]
        x = x[:-7]
        if x.find('_') != -1:
            x = tuple(sorted(x.split('_')))
        else:
            x = tuple(x.split("_"))
        finished[i] = x

    finished = set()
    unfinished = set()

    # count = set()
    #
    # for key, value in samples.items():
    #     count.add(tuple(sorted(key)))
    #     if len(value) != 0:
    #         count.add(tuple(sorted(value)))
    #
    # print len(count)
    # count = 0

    for key, value in samples.items():
        # print key, value
        if 0 <= count < 2000:
            if len(value) >0:
                cur_key = tuple(sorted(key))
                if cur_key not in finished:
                    if not os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/'+'_'.join(cur_key) + '.bowtie'):
                        ENC_multiple_bowtie(cur_key)
                        print '/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/'+'_'.join(cur_key) + '.bowtie'
                        # unfinished.add(cur_key)
                    pass

                    # for cur_k in cur_key:
                    #     if not os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/'+cur_k + '.bowtie'):
                    #         ENC_multiple_bowtie([cur_k])


                cur_value = tuple(sorted(value))
                if cur_value not in finished:
                    if not os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/' + '_'.join(list(cur_value)) + '.bowtie'):
                        print '/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/' + '_'.join(list(cur_value)) + '.bowtie'
                        ENC_multiple_bowtie(cur_value)
                        unfinished.add(cur_value)
                    pass

                    # for cur_v in cur_value:
                    #     if not os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/'+cur_v + '.bowtie'):
                    #         ENC_multiple_bowtie([cur_v])

            else:
                cur_key = tuple(sorted(key))
                if cur_key not in finished:
                    if not os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/' + '_'.join(cur_key) + '.bowtie'):
                        print '/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/' + '_'.join(cur_key) + '.bowtie'
                        ENC_multiple_bowtie(cur_key)
                        unfinished.add(cur_key)
                        pass
                    pass

                    # for cur_k in cur_key:
                    #     if not os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/'+cur_k + '.bowtie'):
                    #         ENC_multiple_bowtie([cur_k])

        if len(value) > 0:

            count += 1
        else:
            count += 1
        results.append(('_'.join(sorted(key)), '_'.join(sorted(value))))
    print count
    result_df = pd.DataFrame(results)

    result_df.columns = ['sample', 'input']

    result_df.to_csv('ENC_H3K4me3_sample_pairs.csv', index=None)

    print len(finished)
    #
    print len(unfinished)

    return results

    # step 3 run bowtie

    # step 4 run danpos




# GEO_sample_input_pair("H3K4me3_GEO_search.csv", "H3K4me3_GEO_metadata.txt")

# files = [x for x in os.listdir("./") if x.endswith('.fastq')]
#
# for f in files:
#     f = f.replace('.fastq', '')
#     SRR_single(f, 1, "Mus musculus")

# pair = open('pairSRR.txt', 'r')
# pair_info = [x.strip() for x in pair.readlines()]
# pair.close()
#
# single = open('singleSRR.txt', 'r')
# single_info =[x.strip() for x in single.readlines()]
# single.close()
#
# # for id in pair_info:
# #     SRR_pair(id,0)
#     # break
#
# for id in single_info[100:400]:
#     SRR_single(id,0)
#     # break










samples = bowtie_ENC('H3K4me3_ENC_metadata.tsv')

print len(samples)


# File accession
# ENCFF240OHP     ENCLB263HJZ                        2                    2_2
# ENCFF723UTH     ENCLB333UZF                        2                    2_1
# ENCFF521WIC     ENCLB263HJZ                        2                    2_2
# set(['ENCFF671ROB', 'ENCFF769IGM', 'ENCFF000CTV'])

# File accession
# ENCFF671ROB     ENCLB399AFA                        3                    3_2
# ENCFF769IGM     ENCLB399AFA                        3                    3_2
# ENCFF000CTV     ENCLB695AMT                        3                    3_1
# set(['ENCFF000CDX', 'ENCFF399MSH'])