## This is module is used to check running result from the standard out and standard error from different analysis software.


import os, pandas as pd, numpy as np, json, urllib2
from collections import defaultdict

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

def bowtie_results(path, surffix='.e', pairs='ENC_H3K4me3_sample_pairs.csv'):
    """
    :param path: the path containing bowtie results in standard error
           surffix: usurally is named as '.error', if not please specify
    :return: a dataframe table that containing error reports, and will print out the list if the error report is empty
    which happened a lot when using sra toolkit
    """
    pairs_df = pd.read_csv(pairs, index_col=0)

    valid_bowties = set()
    for sample in pairs_df.index:
        valid_bowties.add(tuple(sample.split('_')))
        if pd.isnull(pairs_df.ix[sample, 'input']) or len(pairs_df.ix[sample, 'input']) == 0:
            continue
        else:
            valid_bowties.add(tuple(pairs_df.ix[sample, 'input'].split('_')))

    error_files = [x for x in os.listdir(path) if x.find(surffix) != -1]

    correct_error_files = set()
    correct_error_tuple = set()

    for i in range(len(error_files)):
        x = error_files[i]
        x = x[:x.find('.')]
        if x.find('_') != -1:
            x = tuple(sorted(x.split('_')))
        else:
            x = tuple(x.split("_"))
        if x in valid_bowties:
            correct_error_files.add(error_files[i])
            correct_error_tuple.add(x)

    print len(correct_error_files)
    # print len(valid_bowties)
    # for b in valid_bowties:
    #     if b not in correct_bowtie:
    #         print b

    finished = [x for x in os.listdir('/archive2/tmhbxx3/H3K4me3/ENCODE_sample_with_input/bowtie2/') if
                x.endswith('bowtie')]
    for i in range(len(finished)):
        x = finished[i]
        x = x[:-7]
        if x.find('_') != -1:
            x = tuple(sorted(x.split('_')))
        else:
            x = tuple(x.split("_"))

        if x not in valid_bowties:
            # os.system("mv ../bowtie/"+finished[i]+" ../bowtie_old")
            print x
    # finished = set(finished)


    # path = path + '/' if not path.endswith("/") else path
    #
    results = {}
    failed = []

    errors = set()

    for error_file in correct_error_files:
        f = open(path+error_file, "r")
        info = f.readlines()
        f.close()
        error = True
        filename = error_file[:error_file.find(surffix)]
        for line in info:
            if line.startswith("# reads with at least one reported alignment:"):
                error = False
                percentage = line[line.find("(")+1:line.find(")")]
            if line.startswith("# reads processed:"):
                reads_number = line.split()[-1]
            if line.startswith("# reads with alignments suppressed due to -m:"):
                ununique_reads = line[line.find("(")+1:line.find(")")]
            if line.find("storage exhausted while writing file within file system module") != -1:
                errors.add(error_file)
            elif line.find("fastq-dump.2.8.0")!= -1:
                errors.add(error_file)

        # if filename == 'ENCFF000VJK':
        #     print error

        if filename in pairs_df.index:
            type = 'sample'
        elif filename in pairs_df['input'].values:
            type = 'input'

        if error:
            failed.append(filename)
        else:
            if filename not in results:
                results[filename] = (filename, percentage, ununique_reads, reads_number, type)
            else:
                if float(results[filename][1][:-1]) < float(percentage[:-1]):
                    results[filename] = (filename, percentage, ununique_reads, reads_number, type)
    final_results = []
    for key,value in results.items():
        final_results.append(value)

    df = pd.DataFrame(final_results, columns=['sample_ID','percentage', 'ununique_percentage', 'reads_number', 'type'])

    re_failed = []
    for fail in failed:
        if fail not in df['sample_ID'].unique():
            re_failed.append(fail)
    print re_failed
    print errors, len(errors)
    df.to_csv('bowtie_results.csv', index=False)
    return df

def bowtie_results_GEO(path, surffix='.e'):
    """
    :param path: the path containing bowtie results in standard error
           surffix: usurally is named as '.error', if not please specify
    :return: a dataframe table that containing error reports, and will print out the list if the error report is empty
    which happened a lot when using sra toolkit
    """
    GSM_SRR_map = defaultdict(set)
    SRR_GSM_map = {}

    srr_df = pd.read_csv('GEO_sample_input_meta.txt', sep='\t')

    for i in range(srr_df.shape[0]):
        GSM = srr_df.ix[i, 'SampleName']
        GSM2 = srr_df.ix[i, 'GSM_ID']
        if not GSM.startswith('GSM'):
            GSM = None
        if pd.isnull(GSM2):
            GSM2 = None
        if (GSM is not None and GSM2 is not None) and GSM != GSM2:
            print GSM, GSM2, srr_df.ix[i, 'Run_ID']
            continue
        if GSM2 is not None:
            final_gsm = GSM2
        elif GSM is not None:
            final_gsm = GSM
        GSM_SRR_map[final_gsm.strip()].add(srr_df.ix[i, 'Run_ID'])
        SRR_GSM_map[srr_df.ix[i, 'Run_ID']] = final_gsm.strip()

    gsm_df = pd.read_csv('GEO_sample_input.csv')
    pairs = []
    for i in range(gsm_df.shape[0]):
        sample = gsm_df.ix[i, 'Data_ID']
        input = gsm_df.ix[i, 'Input']

        samples = GSM_SRR_map[sample.strip()]
        if len(samples) > 1:
            name = ''
            for s in samples:
                name = name + s + '_'
            name = name[:-1]
            name = '_'.join(sorted(name.split('_')))

            samples = set()
            samples.add(name)

        if len(samples) == 0:
            print gsm_df.ix[i, 'Data_ID'], GSM_SRR_map[gsm_df.ix[i, 'Data_ID'].strip()]
            continue

        inputs = set()
        for f in [x.strip() for x in input.split(',')]:
            if len(GSM_SRR_map[f]) > 1:
                name = ''
                for j in GSM_SRR_map[f]:
                    name = name + j + '_'
                name = name[:-1]
                name = '_'.join(sorted(name.split('_')))
                inputs.add(name)
            else:
                inputs = inputs.union(GSM_SRR_map[f])

        pairs.append((gsm_df.ix[i, 'Data_ID'].strip(), tuple(samples), tuple(inputs)))
    #

    results = {}
    error_files = [x for x in os.listdir(path) if x.find(surffix) != -1]
    error_file_map = {}
    errors =set()
    failed = []
    for error_file in error_files:
        filename = error_file[7:error_file.find(surffix)]
        # print filename
        error_file_map[filename] = error_file
    for p in pairs:
        sample_id, samples, inputs = p

        total_reads = 0
        total_unique = 0
        total_ununique = 0
        for sample in samples:
            ss = sample.split('_')
            for s in ss:
                if s not in error_file_map:
                    print s
                    continue

                f = open(path + error_file_map[s], "r")
                info = f.readlines()
                f.close()
                error = True

                for line in info:
                    if line.startswith("# reads with at least one reported alignment:"):
                        error = False
                        percentage = float(line[line.find("(") + 1:line.find(")")-1])/100
                    if line.startswith("# reads processed:"):
                        reads_number = int(line.split()[-1])
                    if line.startswith("# reads with alignments suppressed due to -m:"):
                        ununique_reads = float(line[line.find("(") + 1:line.find(")")-1])/100
                    if line.find("storage exhausted while writing file within file system module") != -1:
                        errors.add(error_file_map[s])
                    elif line.find("fastq-dump.2.8.0") != -1:
                        errors.add(error_file_map[s])

                total_reads += reads_number
                total_unique += int(reads_number * percentage)
                total_ununique += int(reads_number * ununique_reads)

                if error:
                    failed.append(error_file_map[s])
                    break
        if not error:
            if sample_id not in results:
                results[sample_id] = (sample_id, total_reads, total_unique, total_ununique)
            else:
                if results[sample_id][2] < total_unique:
                    results[sample_id] = (sample_id, total_reads, total_unique, total_ununique)

    final_results = []
    for key,value in results.items():
        final_results.append(value)

    df = pd.DataFrame(final_results, columns=['sample_ID','total_reads', 'unique_reads', 'ununique_reads'])

    re_failed = []
    for fail in failed:
        if fail not in df['sample_ID'].unique():
            re_failed.append(fail)
    print re_failed
    print errors, len(errors)

    error_out = open('srr_errors.txt', 'w')
    for e in errors:
        error_out.write(e[e.find('SRR'):e.find('.')] + '\n')
    error_out.close()

    df.to_csv('bowtie_GEO_results.csv', index=False)
    return df


def seq_depth(meta_table, bowtie_table):
    meta_df = pd.read_csv(meta_table, index_col=0)
    bowtie_df = pd.read_csv(bowtie_table, index_col=0)
    reads_length = []
    depth = []

    for index in bowtie_df.index:
        if index.find("_") == -1:
            if index in meta_df:
                reads_length.append(meta_df.ix[index, 'Read length'])
                depth.append(int(bowtie_df.ix[index, 'reads_number'])*int(meta_df.ix[index, 'Read length'])*1.0/3095677412)
            elif ENC_biosample_meta(index):
                cur_json = ENC_biosample_meta(index)
                reads_length.append(cur_json['read_length'])
                depth.append(int(bowtie_df.ix[index, 'reads_number']) * int(
                    cur_json['read_length']) * 1.0 / 3095677412)
        else:
            indexes = index.split("_")
            total = 0
            cur_reads_length = []
            for j in indexes:
                # print (j in meta_df.index)
                if ENC_biosample_meta(j):
                    cur_json = ENC_biosample_meta(j)
                    if 'read_length' in cur_json:
                        cur_reads_length.append(cur_json['read_length'])
                        total += int(cur_json['read_count']) * int(cur_json['read_length'])
                    else:
                        print j, cur_json.keys()
            print indexes
            cur_reads_length = [str(x) for x in cur_reads_length]
            reads_length.append("_".join(cur_reads_length))
            depth.append(total*1.0/3095677412)

    bowtie_df['Read length'] = reads_length

    bowtie_df['Depth'] = depth

    bowtie_df.to_csv('bowtie_ENC_results.csv')

bowtie_results('./', surffix='.e')
# seq_depth('all_ENC_metadata.csv', 'bowtie_results.csv')

# J = ENC_biosample_meta('ENCFF000CFR')
# print J['read_length']
#
# [u'read_count', u'controlled_by', u'read_length_units', u'file_type', u'fastq_signature', u'accession', u'dataset', u'quality_metrics', u'href', u'technical_replicates', u'alternate_accessions', u'file_size', u'aliases', u'submitted_file_name', u'dbxrefs', u'uuid', u'file_format', u'read_length', u'content_md5sum', u'schema_version', u'platform', u'flowcell_details', u'replicate', u'@context', u'status', u'run_type', u'biological_replicates', u'submitted_by', u'award', u'lab', u'output_type', u'@id', u'audit', u'output_category', u'superseded_by', u'title', u'date_created', u'md5sum', u'@type']