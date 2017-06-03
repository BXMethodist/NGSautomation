"""
This module is used to get biosample, and biological replicate information from Encode Restful API
"""

import urllib2, json, pandas as pd, numpy as np
from collections import defaultdict

def ENC_experiment_meta(exp_id):
    url = "https://www.encodeproject.org/experiments/"+exp_id+"/?format=json"
    response = urllib2.urlopen(url)
    data = json.loads(response.read())
    return data

def ENC_biosample_meta(biosample_id, type):
    # print biosample_id
    if biosample_id is None or pd.isnull(biosample_id):
        return None

    url = "https://www.encodeproject.org/"+type+"/" + biosample_id + "/?format=json"
    try:
        response = urllib2.urlopen(url)
        data = json.loads(response.read())
    except:
        print url
        return None
    return data

def add_library_biosample_metadata(table=None):
    if table is None:
        url = 'https://www.encodeproject.org/metadata/type=Experiment&files.file_type=fastq/metadata.tsv'
        df = pd.read_csv(url, sep='\t', dtype=str, index_col=0)
    else:
        df = pd.read_csv(table, index_col=0)

    df = df[(df['Assay'] == 'RNA-seq') | (df['Assay'] == 'ChIP-seq') | (df['Biosample organism'] == 'Homo sapiens')]

    library_ids = []
    biosample_ids = []
    biological_replicate_ids = []
    technical_replicates_ids = []
    treatments = []

    for exp_id in df.index:
        cur_json = ENC_experiment_meta(exp_id)

        try:
            library_id = cur_json['replicate']['library']['accession']
        except:
            print exp_id
            library_id = None
        try:
            biosample_id = cur_json['replicate']['library']['biosample'].replace('biosamples', '').replace('/', '').strip()
        except:
            print 'biosample', exp_id
            biosample_id = None
        try:
            biological_replicate_id = cur_json['biological_replicates'][0]
            # print biological_replicate_id
        except:
            print 'bioreplicate', exp_id
            biological_replicate_id = None
        try:
            technical_replicates_id = cur_json['technical_replicates'][0]
        except:
            print 'tech', exp_id
            technical_replicates_id = None
        if biosample_id is not None:
            cur_treatment = get_treatment(biosample_id, 'biosamples')
        else:
            print 'treatment', exp_id
            cur_treatment = None

        library_ids.append(library_id)
        biosample_ids.append(biosample_id)
        biological_replicate_ids.append(biological_replicate_id)
        technical_replicates_ids.append(technical_replicates_id)
        treatments.append(cur_treatment)

    df['Library ID'] = library_ids
    df['Biosample ID'] = biosample_ids
    df['Biological Replicate ID'] = biological_replicate_ids
    df['Technical Replicate ID'] = technical_replicates_ids
    df['treatment'] = treatments

    df.to_csv('all_ENC_metadata.csv', encoding='utf-8')

    return

def add_donor(table):
    df = pd.read_csv(table, index_col=0)
    donors_ID = {}

    for biosample in df['Biosample ID'].unique():
        if pd.isnull(biosample):
            continue
        cur_json = ENC_biosample_meta(biosample, 'biosamples')
        if cur_json is None:
            donors_ID[biosample] = None
        try:
            donors_ID[biosample] = cur_json['donor']['accession']
        except:
            donors_ID[biosample] = None

    results = []

    for i in range(df.shape[0]):
        if df.ix[i, 'Biosample ID'] not in donors_ID:
            results.append(None)
        else:
            results.append(donors_ID[df.ix[i, 'Biosample ID']])

    df['donor_ID'] = results

    df.to_csv('all_ENC_metadata.csv')

    return df

def add_assay_title(table):
    df = pd.read_csv(table, index_col=0)
    assay_titles = {}

    for exp in df['Experiment accession'].unique():
        if pd.isnull(exp):
            continue
        cur_json = ENC_biosample_meta(exp, 'experiments')
        if cur_json is None:
            assay_titles[exp] = None
        try:
            assay_titles[exp] = cur_json['assay_title']
        except:
            assay_titles[exp] = None
        print assay_titles

    results = []

    for i in range(df.shape[0]):
        if df.ix[i, 'Experiment accession'] not in assay_titles:
            results.append(None)
        else:
            results.append(assay_titles[df.ix[i, 'Experiment accession']])

    df['assay_title'] = results

    df.to_csv('all_ENC_metadata.csv')

    return df

def add_strand_specificity(table):
    df = pd.read_csv(table, index_col=0)
    strand_specificity = defaultdict(dict)
    count = 0
    for exp in df['Experiment accession'].unique():
        count += 1
        if pd.isnull(exp):
            continue
        cur_json = ENC_biosample_meta(exp, 'experiments')
        if cur_json is None:
            strand_specificity[exp] = None
        else:
            for key in range(len(cur_json['replicates'])):
                cur_rep = cur_json['replicates'][key]
                # print cur_rep['biological_replicate_number'], key
                strand_specificity[exp][cur_rep['biological_replicate_number']] = cur_json['replicates'][key]['library']['strand_specificity']
        # if count > 50:
        #     break
    # print strand_specificity
    # print type(strand_specificity[strand_specificity.keys()[0]][strand_specificity[strand_specificity.keys()[0]].keys()[0]])
    results = []

    for i in range(df.shape[0]):
        if df.ix[i, 'Experiment accession'] not in strand_specificity:
            # print 'no strand specificity!'
            results.append(None)
        else:
            # print df.ix[i, 'Biological Replicate ID']
            # print type(df.ix[i, 'Biological Replicate ID'])
            if df.ix[i, 'Biological Replicate ID'] in strand_specificity[df.ix[i, 'Experiment accession']]:
                results.append(strand_specificity[df.ix[i, 'Experiment accession']][df.ix[i, 'Biological Replicate ID']])
                # print 'found strand_specificity'
            else:
                results.append(None)
                print df.ix[i, 'Biological Replicate ID'], df.ix[i, 'Experiment accession']

    df['strand_specificity'] = results

    df.to_csv('all_ENC_metadata.csv', encoding='utf-8')

    return df

def get_treatment(sample1, type):
    # print sample1
    a = ENC_biosample_meta(sample1, type)
    if a is None:
        return None

    # if 'description' in a:
    #     return a['description']
    if 'summary' in a:
        return a['summary']
    else:
        return None

def same_treatments(sample1, sample2):
    a = ENC_biosample_meta(sample1, 'biosamples')
    b = ENC_biosample_meta(sample2, 'biosamples')
    if a is None or b is None:
        return False
    return a['summary'] == b['summary']

def pandas_compare(feature1, feature2):
    return feature1 == feature2 or (pd.isnull(feature1) and pd.isnull(feature2))

def pair_ENC(table, feature1, type1, type2):
    """
    :param table: organized ENC metadata with biosample ID, donor ID, type1 NGS, type2 NGS
    :param feature1:
    :param type1:
    :param feature2:
    :param type2:
    :return:
    """
    # samples, groups, samples_to_groups = ENC_samples(table)

    df = pd.read_csv(table, index_col=0)

    if feature1 is None or len(feature1) ==0:
        df1 = df[df['assay_title'] == type1].copy()
    else:
        df1 = df[(df['assay_title'] == type1) & (df['Experiment target'] == feature1)].copy()

    # get the candidates of target sample's based on donor ID and
    donor_treatment = defaultdict(set)
    donor_sample_treatment = set()
    # print len(df1['Experiment accession'].unique())
    for exp in df1['Experiment accession'].unique():
        replicates = df[df['Experiment accession'] == exp]['Biological replicate(s)'].unique()
        for rep in replicates:
            cur_rows = df1[(df1['Experiment accession'] == exp) & (df1['Biological replicate(s)'] == rep)]
            cur_donors = cur_rows['donor_ID'].unique()
            # print cur_donors, type(cur_donors)
            if pd.isnull(cur_donors) or len(list(cur_donors)) == 0:
                print exp
            else:
                cur_donor = cur_donors[0]
                treatments = set()
                for biosample_id in cur_rows['Biosample ID'].unique():
                    cur_treatment = get_treatment(biosample_id, 'biosamples')
                    if cur_treatment is not None:
                        treatments.add(cur_treatment)
                # print treatments
                treatments = ';'.join(sorted(list(treatments)))
                # print treatments

                donor_treatment[(cur_donor, treatments)].add((exp, str(rep)))

                donor_sample_treatment.add((cur_donor, tuple(cur_rows['Biosample ID'].unique()), treatments, 'first'))

    # target_donor_biosample_treatment = set()
    target_donor_treatments = {}
    for type in type2:
        df2 = df[(df['assay_title'] == type) & (df['Biosample organism'] == 'Homo sapiens')].copy()

        print len(donor_treatment)
        target_donor_treatment = defaultdict(set)
        for exp in df2['Experiment accession'].unique():
            replicates = df[df['Experiment accession'] == exp]['Biological replicate(s)'].unique()
            for rep in replicates:
                cur_rows = df2[(df2['Experiment accession'] == exp) & (df2['Biological replicate(s)'] == rep)]
                cur_donors = cur_rows['donor_ID'].unique()
                # print cur_donors, type(cur_donors)
                if pd.isnull(cur_donors) or len(list(cur_donors)) == 0:
                    print exp
                else:
                    cur_donor = cur_donors[0]
                    print cur_donor, 'target'
                    treatments = set()
                    for biosample_id in cur_rows['Biosample ID'].unique():
                        cur_treatment = get_treatment(biosample_id, 'biosamples')
                        print cur_treatment, 'target'
                        if cur_treatment is not None:
                            treatments.add(cur_treatment)
                    # print treatments
                    treatments = ';'.join(sorted(list(treatments)))
                    # print treatments

                    target_donor_treatment[(cur_donor, treatments)].add((exp, str(rep)))
                    donor_sample_treatment.add(
                        (cur_donor, tuple(cur_rows['Biosample ID'].unique()), treatments, 'second'))
        target_donor_treatments[type] = target_donor_treatment
    results = []

    for key, value in donor_treatment.items():
        cur_result = [key[0], key[1]]
        for type, target_donor_treatment in target_donor_treatments.items():
            if key in target_donor_treatment:
                if ';'.join(sorted(['_'.join(x) for x in value])) not in cur_result:
                    samples1 = ';'.join(sorted(['_'.join(x) for x in value]))
                    cur_result.append(samples1)
        if ';'.join(sorted(['_'.join(x) for x in value])) not in cur_result:
            continue
        for type, target_donor_treatment in target_donor_treatments.items():
            if key in target_donor_treatment:
                samples2 = ';'.join(sorted(['_'.join(x) for x in target_donor_treatment[key]]))
                cur_result.append(samples2)
            else:
                cur_result.append('')

        results.append(cur_result)

    # bio_conditions = set()
    #
    # for key in donor_treatment.keys():
    #     bio_conditions.add(key)
    # for key in target_donor_treatment.keys():
    #     bio_conditions.add(key)

    # bio_df = pd.DataFrame(list(donor_sample_treatment))
    # bio_df.columns=['donor_ID', "biosample_id" 'treatments', 'feature']
    # bio_df.to_csv('donors_treatments.csv', index=None, encoding='utf-8')

    result_df = pd.DataFrame(results)
    result_df.columns = ['donor_ID', 'treatments', type1+'_samples'] + [type+'_samples' for type in target_donor_treatments.keys()]

    filename = feature1+'_'+type1+"_"+'_'.join(type2)+'.csv'
    filename = filename.replace(' ', "_")

    result_df.to_csv(filename, index=None, encoding='utf-8')
    return

    pairs = defaultdict(set)

    count1 = 0
    count2 = 0
    count3 = 0
    count4 = 0









    for i in range(df1.shape[0]):
        cur_same = []
        cur_same_rep = []
        cur_same_sample = []
        cur_same_sample_rep = []

        candidates = df[(df['donor_ID'] == df1.ix[i, 'donor_ID']) & (df['assay_title'].isin(type2))]
        for j in range(candidates.shape[0]):
            if same_treatments(candidates.ix[j, 'Biosample ID'], df1.ix[i, 'Biosample ID']) and \
                    pandas_compare(candidates.ix[j, 'Biosample term id'], df1.ix[i, 'Biosample term id']) and \
                    pandas_compare(candidates.ix[j, 'Biosample term name'], df1.ix[i, 'Biosample term name']) and \
                    pandas_compare(candidates.ix[j, 'Biosample type'], df1.ix[i, 'Biosample type']) and \
                    pandas_compare(candidates.ix[j, 'Biosample life stage'], df1.ix[i, 'Biosample life stage']) and \
                    pandas_compare(candidates.ix[j, 'Biosample sex'], df1.ix[i, 'Biosample sex']) and \
                    pandas_compare(candidates.ix[j, 'Biosample Age'], df1.ix[i, 'Biosample Age']) and \
                    pandas_compare(candidates.ix[j, 'Biosample organism'], df1.ix[i, 'Biosample organism']):
                if pandas_compare(candidates.ix[j, 'Biosample ID'], df1.ix[i, 'Biosample ID']) and \
                        pandas_compare(candidates.ix[j, 'Biological Replicate ID'],df1.ix[i, 'Biological Replicate ID']):
                    cur_same_sample_rep.append(candidates.iloc[j].name)
                elif pandas_compare(candidates.ix[j, 'Biosample ID'], df1.ix[i, 'Biosample ID']):
                    cur_same_sample.append(candidates.iloc[j].name)
                elif pandas_compare(candidates.ix[j, 'Biological Replicate ID'], df1.ix[i, 'Biological Replicate ID']):
                    cur_same_rep.append(candidates.iloc[j].name)
                else:
                    cur_same.append(candidates.iloc[j].name)

        added = False

        if len(cur_same_sample_rep) > 0:
            for sample in cur_same_sample_rep:
                if df1.iloc[i].name.strip() in samples_to_groups and sample.strip() in samples_to_groups:
                    if samples_to_groups[sample.strip()] not in pairs[samples_to_groups[df1.iloc[i].name.strip()]]:
                        pairs[samples_to_groups[df1.iloc[i].name.strip()]].add(samples_to_groups[sample.strip()])
                        added = True
                        count1 += 1

        if len(cur_same_sample) > 0:
            if not added:
                for sample in cur_same_sample:
                    if df1.iloc[i].name.strip() in samples_to_groups and sample.strip() in samples_to_groups:
                        if samples_to_groups[sample.strip()] not in pairs[samples_to_groups[df1.iloc[i].name.strip()]]:
                            pairs[samples_to_groups[df1.iloc[i].name.strip()]].add(samples_to_groups[sample.strip()])
                            added = True
                            count2 += 1

        if len(cur_same_rep) > 0:
            if not added:
                for sample in cur_same_rep:
                    if df1.iloc[i].name.strip() in samples_to_groups and sample.strip() in samples_to_groups:
                        if samples_to_groups[sample.strip()] not in pairs[samples_to_groups[df1.iloc[i].name.strip()]]:
                            pairs[samples_to_groups[df1.iloc[i].name.strip()]].add(samples_to_groups[sample.strip()])
                            count3 += 1
                            added = True


        if len(cur_same) > 0:
            if not added:
                for sample in cur_same:
                    if df1.iloc[i].name.strip() in samples_to_groups and sample.strip() in samples_to_groups:
                        if samples_to_groups[sample.strip()] not in pairs[samples_to_groups[df1.iloc[i].name.strip()]]:
                            pairs[samples_to_groups[df1.iloc[i].name.strip()]].add(samples_to_groups[sample.strip()])
                            added = True
                            count4 += 1
        # else:
        #     for sample in cur_same_sample_rep:
        #         pairs.add((samples_to_groups[df1.iloc[i].name.strip()], samples_to_groups[sample.strip()]))

    # df1['Paired_RNA_seq'] = results

    # df1.to_csv('H3K4me3_ENC_paired_chip_rna_seq.csv')
    print pairs
    print len(pairs)

    total = 0

    output = open('H3K4me3_chipseq_RNA_seq_pair_ENC.csv', 'w')

    for key, value in pairs.items():
        total += len(value)
        output.write(key+','+';'.join(list(value)) + '\n')
    print total
    print count1, count2, count3, count4

    output.close()
    return

def ENC_samples(table):
    df = pd.read_csv(table, index_col=0)
    groups = {}
    samples = {}
    samples_to_groups = {}

    for experiment in df['Experiment accession'].unique():
        replicates = df[df['Experiment accession'] == experiment]['Biological replicate(s)'].unique()
        for rep in replicates:
            # try:
            cur_samples = df[(df['Experiment accession'] == experiment) & (df['Biological replicate(s)'] == rep)].index
            cur_samples = set(cur_samples)

            valid = True

            known_samples = cur_samples.copy()

            for id in known_samples:
                if not pd.isnull(df.ix[id, 'Paired with']):
                    cur_samples.add(df.ix[id, 'Paired with'].strip())

            cur_inputs = set()

            for id in cur_samples:
                if id in df.index and not pd.isnull(df.ix[id, 'Controlled by']):
                    cur_inputs = [input.replace('/', '').replace('files', '').strip() for input in
                                  df.ix[id, 'Controlled by'].split(',')]
                else:
                    cur_inputs = []

            for sample in cur_samples:
                # print df.loc[sample].name
                if sample not in df.index:
                    print sample, 'is not here!'
                    samples[sample.strip()] = None
                    valid = False
                else:
                    samples[sample.strip()] = ENC_file(df.loc[sample].name,
                                                       df.ix[sample, 'Experiment accession'],
                                                       df.ix[sample, 'Assay'],
                                                       df.ix[sample, 'assay_title'],
                                                       df.ix[sample, 'Library ID'],
                                                       df.ix[sample, 'Biosample ID'],
                                                       df.ix[sample, 'Biosample term id'],
                                                       df.ix[sample, 'Biosample term name'],
                                                       df.ix[sample, 'Biosample type'],
                                                       df.ix[sample, 'Biosample life stage'],
                                                       df.ix[sample, 'Biosample sex'],
                                                       df.ix[sample, 'Biosample Age'],
                                                       df.ix[sample, 'Biosample organism'],
                                                       df.ix[sample, 'Experiment target'],
                                                       df.ix[sample, 'Read length'],
                                                       df.ix[sample, 'Run type'],
                                                       df.ix[sample, 'Paired with'],
                                                       df.ix[sample, 'Paired end'],
                                                       df.ix[sample, 'Controlled by'],
                                                       df.ix[sample, 'donor_ID'],
                                                       df.ix[sample, 'Biological Replicate ID'])
            for sample in cur_inputs:
                # print sample
                if sample not in df.index:
                    print sample, 'is not here!'
                    samples[sample.strip()] = None
                    valid = False
                else:
                    samples[sample.strip()] = ENC_file(df.loc[sample].name,
                                                       df.ix[sample, 'Experiment accession'],
                                                       df.ix[sample, 'Assay'],
                                                       df.ix[sample, 'assay_title'],
                                                       df.ix[sample, 'Library ID'],
                                                       df.ix[sample, 'Biosample ID'],
                                                       df.ix[sample, 'Biosample term id'],
                                                       df.ix[sample, 'Biosample term name'],
                                                       df.ix[sample, 'Biosample type'],
                                                       df.ix[sample, 'Biosample life stage'],
                                                       df.ix[sample, 'Biosample sex'],
                                                       df.ix[sample, 'Biosample Age'],
                                                       df.ix[sample, 'Biosample organism'],
                                                       df.ix[sample, 'Experiment target'],
                                                       df.ix[sample, 'Read length'],
                                                       df.ix[sample, 'Run type'],
                                                       df.ix[sample, 'Paired with'],
                                                       df.ix[sample, 'Paired end'],
                                                       df.ix[sample, 'Controlled by'],
                                                       df.ix[sample, 'donor_ID'],
                                                       df.ix[sample, 'Biological Replicate ID'])

            cur_samples_id = "_".join([x.strip() for x in cur_samples])
            cur_inputs_id = "_".join([y.strip() for y in cur_inputs])

            if valid:
                groups[cur_samples_id] = cur_inputs_id

                for sample in cur_samples:
                    samples_to_groups[sample.strip()] = cur_samples_id
                for sample in cur_inputs:
                    samples_to_groups[sample.strip()] = cur_inputs_id
            # except:
            #     continue
    # step 2 cat the fastq in each group and run bowtie
    print len(samples_to_groups)
    # print samples_to_groups.keys()[:100]

    return samples, groups, samples_to_groups

class ENC_sample:
    def __init__(self, sample_name, samples, inputs):
        self.sample_name = sample_name
        self.samples = samples
        self.inputs = inputs

class ENC_file:
    def __init__(self,
                 accession,
                 experiment,
                 assay,
                 assay_title,
                 library_id,
                 biosample_id,
                 biosample_term_id,
                 biosample_term_name,
                 biosample_type,
                 biosample_life_stage,
                 biosample_sex,
                 biosample_age,
                 biosample_organism,
                 experiment_target,
                 read_length,
                 run_type,
                 paired_with,
                 paired_end,
                 controlled_by,
                 donor_ID,
                 replicate_id):
        self.accession = accession
        self.experiment = experiment
        self.assay = assay
        self.assay_title = assay_title
        self.library_id = library_id
        self.biosample_id = biosample_id
        self.biosample_term_id = biosample_term_id
        self.biosample_term_name = biosample_term_name
        self.biosample_type = biosample_type
        self.biosample_life_stage = biosample_life_stage
        self.biosample_sex = biosample_sex
        self.biosample_age = biosample_age
        self.biosample_organism = biosample_organism
        self.experiment_target = experiment_target
        self.read_length = read_length
        self.run_type = run_type
        self.paired_with = paired_with
        self.paired_end = paired_end
        self.controlled_by = controlled_by
        self.donor_ID = donor_ID
        self.replicate_id = replicate_id

#
#
# data = ENC_experiment_meta('ENCSR000DSD')
#
# print data['files'][6]['accession']
#
#
#
# data['replicates'][0]['biological_replicate_number']
# data['replicates'][0]['library']['accession']
# data['replicates'][0]['library']['biosample']['accession']
# data['files'][5]['file_type']

# df = pd.read_csv('all_ENC_metadata.tsv', sep='\t', index_col=0)
# df2 = pd.read_csv('biosample_ids.tsv', sep='\t', index_col=0)
#
# result_df = df.join(df2)
#
# print result_df[result_df['Biological replicate(s)'] != result_df['Biological Replicate ID']]
#
# result_df.to_csv('all_ENC_metadata_biosample.csv')

# add_donor('all_ENC_metadata_biosample.csv')

# add_assay_title('all_ENC_metadata.csv')
add_library_biosample_metadata('all_ENC_metadata.csv')
add_strand_specificity('all_ENC_metadata.csv')
#
# pair_ENC('all_ENC_metadata.csv', '', 'polyA mRNA RNA-seq', ['total RNA-seq'])

# a = ENC_biosample_meta('ENCSR656FIH')
# b = ENC_biosample_meta('ENCBS160AAA')
# print a.keys()
# print a['treatments'] == b['treatments']
#
# print a['replicates'][0]['library']['strand_specificity'] is True

# print same_treatments(None, None)

# add_library_biosample_metadata()
# add_donor('all_ENC_metadata.csv')
# add_assay_title('all_ENC_metadata.csv')

# cur_json = ENC_experiment_meta('ENCFF668QER')
#
# print cur_json['replicate']['library']['biosample'].replace('biosamples', '').replace('/', '').strip()



