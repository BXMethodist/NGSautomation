# This module contains functions to run danpos and related function
import os, pandas as pd
from collections import defaultdict

def danpos_no_input(sample_id, search_df, metadata_df, node_id):
    if sample_id.startswith("GSM"):
        sample_id = metadata_df.ix[sample_id, 'Run_ID']

    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' --smooth_width 0 -c 25000000 --frsz 200 --extend 200 -o ' + os.getcwd() + '/' + sample_id

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

def danpos_ENC_no_input(sample_id):
    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' --smooth_width 0 -c 25000000 --frsz 200 --extend 200 -o ' + os.getcwd() + '/' + sample_id

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

def danpos_input(sample_id, input_id, search_df, metadata_df, node_id):
    if sample_id.startswith("GSM"):
        sample_id = metadata_df.ix[sample_id, 'Run_ID']
        input_id = metadata_df.ix[input_id, 'Run_ID']

    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' --smooth_width 0 -c 25000000 --frsz 200 --extend 200 -o ' + os.getcwd() + '/' + sample_id

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

def danpos_input_ENC(sample_id, input_id, node_id):
    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' --smooth_width 0 -c 25000000 --frsz 200 --extend 200 -o ' + os.getcwd() + '/' + sample_id

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
    # os.system('qsub ' + sample_id + ".pbs")
    return

def danpos_inputs(sample_id):
    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' --smooth_width 0 -c 25000000 --frsz 200 --extend 200 -o ' + os.getcwd() + '/' + sample_id

    cmd = danpos_cmd + sample_id + '.bowtie' +' -b '+ os.getcwd() + '/' + sample_id+'_input' + danpos_parameters

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
    # os.system('qsub ' + sample_id + ".pbs")
    return

def danpos_multi_SRR(sample_id, input_id):
    if isinstance(sample_id, list):
        sample_id = cat_sample_SRR(sample_id)
        pass
    if isinstance(input_id, list):
        input_id = cat_sample_SRR(input_id)
        pass

    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' --smooth_width 0 -c 25000000 --frsz 200 --extend 200 -o ' + os.getcwd() + '/' + sample_id

    cmd = danpos_cmd + sample_id + '.bowtie' + ' -b ' + input_id + ".bowtie" + danpos_parameters

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

def cat_sample_SRR(ids, path='/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/'):
    outputname = '_'.join(ids)
    cmd = 'cat '
    for id in ids:
        cmd += path+id+'.bowtie '
    cmd +='> '+outputname+'.bowtie'
    return outputname

def RunDanpos(search, metadata, sample_input="ENC_H3K4me3_sample_pairs.csv"):
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

    for i in range(sample_input_df.shape[0])[0:100]:
        node_id = nodes[i%6]
        sample_id = sample_input_df.ix[i, 0]
        input_id = sample_input_df.ix[i, 1]

        if input_id.find(';') != -1 or pd.isnull(input_id):
            continue
        if os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/'+sample_id+'.bowtie') \
            and os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/'+input_id+'.bowtie'):
            os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/"+sample_id+'.bowtie '+sample_id+'.bowtie')
            os.system(
                "cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/" + input_id + '.bowtie ' + input_id + '.bowtie')
        else:
            continue

        if pd.isnull(input_id):
            danpos_no_input(sample_id, search_df, metadata_df, node_id)
        else:
            danpos_input(sample_id, input_id, search_df, metadata_df, node_id)

def RunDanposENC(sample_input="ENC_H3K4me3_sample_pairs.csv"):
    """
    :param sample_input: file path for sample input pair
    :return:
    """
    sample_input_df = pd.read_csv(sample_input, index_col=None)

    finished = [x for x in os.listdir('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/') if
                x.endswith('bowtie')]

    ### This is the special part to re-run error bowtie file
    # incomplete_f = open('incomplete.txt', 'r')
    #
    # incomplete =[line.strip() for line in incomplete_f]
    # incomplete_f.close()
    # incomplete_set = set()
    ###

    names = {}
    for i in range(len(finished)):
        x = finished[i]
        x = x[:-7]
        if x.find('_') != -1:
            new_name = '_'.join(tuple(sorted(x.split('_'))))
        else:
            new_name = x
        names[new_name] = x

    # print names

    nodes = [1, 2, 3, 4, 5, 6]

    count = 0
    for i in range(sample_input_df.shape[0]):
        node_id = nodes[i%6]
        sample_id = sample_input_df.ix[i, 'sample']
        input_id = sample_input_df.ix[i, 'input']

    ### This is the special part to re-run error bowtie file
        # if sample_id in incomplete or input_id in incomplete:
        #     incomplete_set.add(sample_id)
        #     incomplete_set.add(input_id)
        #     if sample_id in incomplete:
        #         print sample_id, 'incomplete'
        #     if input_id in incomplete:
        #         print input_id, 'incomplete'
        #     # os.system('rm -r '+sample_id.strip()+'/')
        #     pass
        #     # continue
        # else:
        #     continue
    # print len(incomplete_set)
    ###

    # '''
        if sample_id.find('_') != -1:
            if '_'.join(tuple(sorted(sample_id.split('_')))) in names:
                sample_id = names['_'.join(tuple(sorted(sample_id.split('_'))))]
            else:
                print sample_id
                pass

        if pd.isnull(input_id):
            print sample_id, 'no input'
            ## sample without inputs will not be run
            # if os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/' + sample_id + '.bowtie') and \
            #     (not os.path.isfile('./' + sample_id + '.bowtie')):
            #     os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/" + sample_id + '.bowtie ' + sample_id + '.bowtie')
            #     danpos_ENC_no_input(sample_id)
            # else:
            #     if not os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/' + sample_id + '.bowtie'):
            #         print sample_id, "bowtie not found"
            continue
        else:
            pass

        if input_id.find('_') != -1 and input_id.find(';') == -1:
            if '_'.join(tuple(sorted(input_id.split('_')))) in names:
                input_id = names['_'.join(tuple(sorted(input_id.split('_'))))]
            else:
                print input_id
                pass


        if input_id.find(';') != -1:
            inputs = [x.strip() for x in input_id.split(';')]
            if os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/' + sample_id + '.bowtie') and \
                    (not os.path.isfile('./' + sample_id + '.bowtie')):
                os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/" + sample_id + '.bowtie ' + sample_id + '.bowtie')
                pass
            else:
                print sample_id, 'sample bowtie not found'
                pass
            # os.system('mkdir '+sample_id)
            for input_file in inputs:
                if input_file.find('_') != -1:
                    input_file = names['_'.join(tuple(sorted(input_file.split('_'))))]
                if os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/'+input_file+'.bowtie') and \
                        (not os.path.isfile('./' + input_file + '.bowtie')):
                    os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/" + input_file + '.bowtie ' + '/' + sample_id+'_input/'+ input_file + '.bowtie')
                    pass
                else:
                    print input_file, 'input bowtie not found'

            danpos_inputs(sample_id)
            count += 1

        elif os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/'+sample_id+'.bowtie') \
            and os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/'+input_id+'.bowtie'):

            if (not os.path.isfile('./' + sample_id + '.bowtie')):
                os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/"+sample_id+'.bowtie '+sample_id+'.bowtie')
            if (not os.path.isfile('./' + input_id + '.bowtie')):
                os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/" + input_id + '.bowtie ' + input_id + '.bowtie')
            danpos_input_ENC(sample_id, input_id, node_id)
            count +=1
            pass
        else:
            print sample_id, input_id, 'sample or input bowtie not found'
            continue

    for pbs in [x for x in os.listdir('./') if x.endswith('.pbs')]:
        os.system('qsub '+pbs)
    print count
    # '''

def RunDanpos_multiple_inputs(search, metadata, sample_input="ENC_H3K4me3_sample_pairs.csv"):
    """
    :param sample_input: file path for sample input pair
    :param search: search dataframe from chipseqpair search function
    :param metadata: metadata dataframe from chipseqpair query function
    :return:
    """
    sample_input_df = pd.read_csv(sample_input, header=None, index_col=None)
    search_df = pd.read_csv(search, index_col=None)
    metadata_df = pd.read_csv(metadata, index_col=None, sep='\t')
    ENC_metadata_df = metadata_df.copy()
    ENC_metadata_df = ENC_metadata_df.set_index(['Run_ID'])

    nodes = [1,2,3,4,5, 6]
    paired = {}
    paired_ids = [y.replace('.bowtie', '') for x in
                  os.listdir('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/') if x.find('_') != -1 for y in
                  x.split('_')]
    for pair in [x for x in os.listdir('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/') if x.find('_') != -1]:
        pair = pair.replace('.bowtie','')
        id1, id2 = pair.split('_')
        paired[id1] = id2
        paired[id2] = id1

    finished = set()
    for i in range(sample_input_df.shape[0])[100:]:
        node_id = nodes[i%6]
        sample_id = sample_input_df.ix[i, 0]
        input_id = sample_input_df.ix[i, 1]

        if input_id.find(';') == -1 or pd.isnull(input_id) or sample_id in finished:
            continue
        # if os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/'+sample_id+'.bowtie') \
        #     and os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/'+input_id+'.bowtie'):

            # os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/"+sample_id+'.bowtie '+sample_id+'.bowtie')
            # os.system(
            #     "cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/" + input_id + '.bowtie ' + input_id + '.bowtie')
        if (os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/'+sample_id+'.bowtie')) or \
                (sample_id in paired_ids):

            if sample_id in paired_ids:
                finished.add(sample_id)
                finished.add(paired[sample_id])

                if os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/'+sample_id+'_'+paired[sample_id]+'.bowtie'):
                    sample_id = sample_id +'_' + paired[sample_id]
                else:
                    sample_id = paired[sample_id] + '_' + sample_id
            else:
                finished.add(sample_id)

            os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/" + sample_id + '.bowtie ' + sample_id + '.bowtie')
            os.system('mkdir '+sample_id+'_input')

            input_ids = input_id.split(";")

            finished_input = set()
            for id in input_ids:
                if id in finished_input:
                    continue
                if os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/' + id + '.bowtie') or id in paired_ids:
                    if id in paired_ids:
                        finished_input.add(id)
                        finished_input.add(paired[id])

                        if os.path.isfile('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/' + id + '_' +
                                                paired[id] + '.bowtie'):
                            id = id + '_' + paired[id]
                        else:
                            id = paired[id] + '_' + id
                    else:
                        finished_input.add(id)
                    os.system("cp /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie/" + id + '.bowtie ' + sample_id+'_input/'+id + '.bowtie')
                # else:
                #     print id
                    # try:
                    #     print ENC_metadata_df.ix[id, 'Run type'], 'input'
                    # except:
                    #     print 'no metadata'
            danpos_inputs(sample_id)
                # else:

                #     continue

        else:
            print sample_id
            try:
                print ENC_metadata_df.ix[sample_id, 'Run type'], 'but sample is not processed'
            except:
                print 'no metadata, sample is not there'
            pass



        # if pd.isnull(input_id):
        #     danpos_no_input(sample_id, search_df, metadata_df, node_id)
        # else:
        #     danpos_input(sample_id, input_id, search_df, metadata_df, node_id)

def RunDanpos_GEO_single_input(search, metadata, sample_input="sample_input_pair.csv"):
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

    count = 0

    for i in range(sample_input_df.shape[0])[0:100]:
        node_id = nodes[i%6]
        sample_id_GSM = sample_input_df.ix[i, 0]
        input_id_GSM = sample_input_df.ix[i, 1]

        if input_id_GSM.find(';') != -1 or pd.isnull(input_id_GSM):
            continue

        try:
            sample_id = metadata_df.ix[sample_id_GSM, 'Run_ID']
            if isinstance(sample_id,pd.Series):
                sample_id = sample_id.values

            input_id = metadata_df.ix[input_id_GSM, 'Run_ID']
            if isinstance(input_id, pd.Series):
                input_id = input_id.values
            # print sample_id, input_id
            # print ''
            count += 1

            if isinstance(sample_id, list) or isinstance(input_id, list):
                # pass
                danpos_multi_SRR(sample_id, input_id)
            else:
                if os.path.isfile('/archive2/tmhbxx3/H3K4me3/GEO_sample_with_input/bowtie/' + sample_id + '.bowtie') \
                        and os.path.isfile(
                                            '/archive2/tmhbxx3/H3K4me3/GEO_sample_with_input/bowtie/' + input_id + '.bowtie'):
                    os.system(
                        "cp /archive2/tmhbxx3/H3K4me3/GEO_sample_with_input/bowtie/" + sample_id + '.bowtie ' + sample_id + '.bowtie')
                    os.system(
                        "cp /archive2/tmhbxx3/H3K4me3/GEO_sample_with_input/bowtie/" + input_id + '.bowtie ' + input_id + '.bowtie')
                    # pass
                else:
                    print "bowtie file is missing", sample_id, input_id
                    continue

                if pd.isnull(input_id):
                    danpos_no_input(sample_id, search_df, metadata_df, node_id)
                else:
                    danpos_input(sample_id, input_id, search_df, metadata_df, node_id)

        except:
            print 'metada not found'
            continue


        # print count

def move_danpos_wig(path='./'):
    # os.system('mkdir ../wigs')
    folders = [x for x in os.listdir(path) if os.path.isdir(x) and x.find('input') == -1]
    for folder in folders:
        cur_path = path+folder+'/pooled/'
        os.system('mv '+cur_path+'*.wig /archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/wigs/')

def generate_wig(f):
    name = f[f.rfind('/')+1:]
    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' -u 1 --smooth_width 0 -c 25000000 --frsz 200 --extend 200 ' \
                        '-o ' + os.getcwd() + '/' + name

    cmd = danpos_cmd + f + danpos_parameters

    pbs = open(name + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N danpos_" + name + '\n')
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
    os.system('qsub ' + name + ".pbs")
    return

def RunDanposGEO():
    df = pd.read_csv('GEO_danpos_pair.csv')
    #
    path = '/archive2/tmhbxx3/H3K4me3/GEO_sample_with_input/bowtie/'
    for i in range(df.shape[0])[:300]:
        sample = df.ix[i, 'Run_ID']
        input = df.ix[i, 'Input_Run_ID']

        if not os.path.isfile(path+sample+'.bowtie') or not os.path.isfile(path+input+'.bowtie'):
            continue

        sample_id = df.ix[i, 'Data_ID']
        os.system('mkdir '+sample_id+'_sample')
        os.system('mkdir '+sample_id+'_input')

        os.system('cp ' + path + sample + '.bowtie ' + './' + sample_id + '_sample/')
        os.system('cp ' + path + input + '.bowtie ' + './' + sample_id + '_input/')

        GEO_danpos(sample_id)
    for pbs in [x for x in os.listdir('./') if x.endswith('.pbs')]:
        os.system('qsub '+pbs)

def GEO_danpos(sample_id):
    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' --smooth_width 0 -c 25000000 --frsz 200 --extend 200 -o ' + os.getcwd() + '/' + sample_id

    cmd = danpos_cmd + os.getcwd() + '/' + sample_id + '_sample' + ' -b ' + os.getcwd() + '/' + sample_id + '_input' + danpos_parameters

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
    # os.system('qsub ' + sample_id + ".pbs")
    return



# groups = [x for x in os.listdir('/home/tmhbxx3/archive/Erin/Chip-seq/bowtie/') if x.endswith('.bowtie')]
# for group in groups:
#     path = '/home/tmhbxx3/archive/Erin/Chip-seq/bowtie/'
#     cur_path = path + group
#     generate_wig(cur_path)



# RunDanpos_GEO_single_input("H3K4me3_GEO_search.csv", "H3K4me3_GEO_metadata.txt")
# move_danpos_wig()

# RunDanposENC(sample_input='archived_sample_input_pair.csv')

RunDanposGEO()