import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
### tophat2 directory
### /home/tmhbxx3/tools/tophat-2.1.1/tophat
### /archive/tmhkxc48/ref_data/hg19/bowtie2/hg19

### fastq-dump
### /home/tmhbxx3/tools/sratoolkit/bin/fastq-dump

### hisat2
### /home/tmhbxx3/tools/hisat2-2.0.5/hisat2
### /home/tmhbxx3/archive/ref_data/hg19/hisat2/hg19


def runTopHat2(path, surffix_R1, surffix_R2, file_type, finishedjob=[], tophatIndex=" /archive/tmhkxc48/ref_data/hg19/bowtie2/hg19 "):
    cmd = "tophat2 -p 8 "

    if not path.endswith("/"):
        path += "/"

    tophatpath = path + "tophat_output/"
    fastqpath = path+"FASTQ/"

    fileNames = os.listdir(fastqpath)

    if not os.path.isdir(tophatpath):
        os.system("mkdir "+tophatpath)

    candidates = set(fileNames)

    for name in fileNames:
        if not name.endswith(file_type):
            continue

        if name[:name.find(surffix_R1)] in finishedjob:
            continue

        if name.endswith(surffix_R1):
            reversePair = name[:name.find(surffix_R1)]+surffix_R2
            if reversePair in candidates:
                outputPath = "-o " + tophatpath + name[:name.find(surffix_R1)]
                curCmd = cmd + outputPath + tophatIndex + fastqpath + name + " " + fastqpath + reversePair
                # print curCmd
                os.system(curCmd)
        elif name.endswith(surffix_R2):
            pass
        else:
            outputPath = "-o ./tophat_output/" + name[:name.find(file_type)]
            curCmd = cmd + outputPath + tophatIndex + fastqpath + name
            # print curCmd
            # os.system(curCmd)
    return


def runFastqdump(list, path):
    filenames = os.listdir(path)

    f = open(list, "r")
    ids = [x for x in f.readlines()]
    f.close()
    cmd = "/home/tmhbxx3/tools/sratoolkit/bin/fastq-dump -O /home/tmhbxx3/archive/H3K4me3/GEO_with_input/input/FASTQ --split-3 "
    for name in filenames:
        if name in ids:
            os.system(cmd + path +name)
    return

def moveAndChangeNameForTophat(path):
    '''

    :param path: the directory contains tophat_output
    :return:
    '''

    if not path.endswith("/"):
        path += "/"

    input_path = path + "tophat_output/"
    output_path = path + "bamfiles/"

    if not os.path.isdir(output_path):
        os.system("mkdir "+output_path)

    SRRlist = os.listdir(input_path)

    for SRR_name in SRRlist:
        old_file_location = input_path+SRR_name+"/accepted_hits.bam"
        new_file_location = input_path+SRR_name+"/"+SRR_name+".bam"
        os.system("mv " + old_file_location + " "+ new_file_location)
        os.system("cp "+new_file_location+" "+output_path)
    return

def get_tophatalign_result(path, output, type="paired", SRRinfo="/home/tmhbxx3/scratch/XMLhttp/pickles/", sample_type=None):
    # path is the directory contains all the samples tophat results, it contains a list of folder name with sample's name and inside are tophat output files
    # SRRinfo is the place save SRR object information
    # sample_type is the GSM to sample type csv file.

    results = {}
    samples = os.listdir(path)

    if not path.endswith("/"):
        path += "/"
    for sample in samples:
        file_path = path+sample+"/align_summary.txt"
        sample_obj = open(file_path, "r")
        info = sample_obj.readlines()
        sample_obj.close()

        for line in info:
            print line

        # total_reads = 0

    #     if type=="paired":
    #         for i in range(len(info)):
    #             line = info[i]
    #             line = line.strip()
    #
    #             if line.startswith("Input"):
    #                 line = [x.strip() for x in line.split(":")]
    #                 total_reads = int(line[1])
    #             elif line.startswith("Aligned pairs"):
    #                 line = [x.strip() for x in line.split(":")]
    #                 aligned_reads = int(line[1])
    #
    #                 multiple_align_line = info[i+1]
    #                 multiple_align_info = [x.strip() for x in multiple_align_line.split(":")][1]
    #                 multiple_align = int(multiple_align_info[0])
    #
    #                 discordant_line = info[i+2]
    #                 discordant_info = [x.strip() for x in discordant_line.split()]
    #                 discordant = int(discordant_info[0])
    #
    #                 percentage = (aligned_reads-multiple_align-discordant)*100.0/total_reads
    #
    #                 # print total_reads, percentage
    #                 results[sample] = percentage
    #                 break
    # df = pd.DataFrame.from_dict(results, orient="index")
    # df.columns = ['mapping_percentage']
    #
    # if sample_type is not None:
    #     import pickle
    #     gsmsrr_map_obj = open(SRRinfo+"GSMSRR_map.pkl", "rb")
    #     gsmsrr_map = pickle.load(gsmsrr_map_obj)
    #     gsmsrr_map_obj.close()
    #
    #     srrtogsm = {}
    #     import csv
    #     reader = csv.reader(open(sample_type, 'r'))
    #     srrtotype = {}
    #     for row in reader:
    #         k, v = row
    #         SRRids = gsmsrr_map[k]
    #         for SRRid in SRRids:
    #             srrtogsm[SRRid] = k
    #             srrtotype[SRRid] = v
    #
    #     df["GSM"] = pd.Series(srrtogsm)
    #     df["sample_type"] = pd.Series(srrtotype)
    #
    # print df
    #
    # df.to_csv(output, sep="\t")

    # return df


def danposRecall(k, cutoff):
    wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"

    wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]


    n = 0
    for wig in wigFiles:
        if 180*k <=n< 180*(k+1):
            cmd = "python /archive/tmhkxc48/tools/danposTemp/danpos.py dpeak "+wigPath+wig+" -q "+str(cutoff)+" -f 0 -z 0 -o /home/tmhbxx3/archive/KFH3K4me3/"+str(cutoff)+"cutoff"
            os.system(cmd)
        n+=1

def GeneNameConverter(gtf, gene_name):
    # this is used for replace the gene names in gtf to common gene name
    gene_name_file = open(gene_name, "r")
    gene_name_map = {}
    for line in gene_name_file.readlines():
        geneid, name = line.rstrip().split("\t")
        gene_name_map[geneid] = name
    gene_name_file.close()

    gtf_file = open(gtf, "r")
    out_put_gtf = open(gtf[:-4]+"_common_name.gtf", "w")
    for line in gtf_file.readlines():
        info = line.split("\t")
        gene_name_info = info[-1]
        input_name = gene_name_info[gene_name_info.find("\"")+1:gene_name_info.find(";")-1]

        if input_name in gene_name_map:
            new_name = gene_name_info[:gene_name_info.find("\"")+1] + gene_name_map[input_name] + \
                       gene_name_info[gene_name_info.find(";") - 1:]
            info[-1] = new_name
            new_line ="\t".join(info)
            out_put_gtf.write(new_line)
        else:
            out_put_gtf.write(line)
    out_put_gtf.close()


def get_danpos_profile_gene_list(path="/home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/", geneinfo="/home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.xls"):
    geneinfolist = pd.read_csv(geneinfo,sep="\t")

    geneinfolist = geneinfolist.set_index('mm9.kgXref.geneSymbol', drop=False)

    file_names = [x for x in os.listdir(path) if x.endswith('.txt')]

    for name in file_names:
        genelist_obj = open(path+name, "r")
        gene_names = []
        for line in genelist_obj.readlines():
            line = line.rstrip().capitalize()
            gene_names.append(line)

        genelist_obj.close()
        genelist_info = geneinfolist.ix[gene_names, :]
        genelist_info.to_csv(name+"_profilelist.txt", sep="\t", index=False, float_format='%d')


def cuffdiff_table_selector(table, name, **kwargs):
    ## **kwargs pair of column name and cutoffs
    ##TO DO: create table selecter method for heatmaps and downstream analysis
    df = pd.read_csv(table, sep="\t")
    df['gene_id'] = df['gene_id'].str.lower()
    df = df.set_index(['gene_id'])
    output_names = [name]
    result_df = None

    for key, value in kwargs.items():
        if key != "gene_id":
            output_names.append(key)
            output_names.append(str(value))
        if result_df is None:
            if key == "gene_id":
                # print value
                result_df = df.ix[value, :]
                # print result_df
            if key == "q_value":
                result_df = df[df['q_value'] < value]
            if key == "p_value":
                result_df = df[df['p_value'] < value]
            if key == "log2(fold_change)":
                if value > 0:
                    result_df = df[df['log2(fold_change)'] > value]
                if value < 0:
                    result_df = df[df['log2(fold_change)'] < value]
            if key == "sample_1":
                result_df = df[df['sample_1'] == value]
            if key == "sample_2":
                result_df = df[df['sample_2'] == value]
            if key == "columns":
                result_df = df[value]
        else:
            if key == "gene_id":
                result_df = result_df.ix[value, :]
            if key == "q_value":
                result_df = result_df[result_df['q_value'] < value]
            if key == "p_value":
                result_df = result_df[result_df['p_value'] < value]
            if key == "log2(fold_change)":
                if value > 0:
                    result_df = result_df[result_df['log2(fold_change)'] > value]
                if value < 0:
                    result_df = result_df[result_df['log2(fold_change)'] < value]
            if key == "sample_1":
                result_df = result_df[result_df['sample_1'] == value]
            if key == "sample_2":
                result_df = result_df[result_df['sample_2'] == value]
            if key == "columns":
                result_df = result_df[value]

    output_filename = "_".join(output_names) + ".tsv"

    # del result_df['test_id']
    # del result_df['gene']
    # del result_df['locus']
    # del result_df['test_stat']

    result_df.to_csv(output_filename, sep="\t")

    print output_filename

    return result_df, "_".join(output_names)


def get_exp_cuffdiff(genelist, **kwargs):
    ### **kwargs are pair of groupname and cuffdiff table
    ## to do pd.concat based on row id
    pass


def get_gene_list(table, name, **kwargs):
    df, output_names = cuffdiff_table_selector(table, name, **kwargs)
    output_names = "_".join(output_names) + "_genelist.csv"

    gene_list = df.index.values

    with open(output_names, "w") as output_file:
        for gene in gene_list:
            output_file.write(gene + "\n")


def gene_exp_barplot(fpkm_table, gene_list):
    #TODO !
    notch_list_obj = open(gene_list, "r")

    notch_gene_list = [x.strip().lower() for x in notch_list_obj.readlines()]
    notch_list_obj.close()

    chars = {"gene_id": notch_gene_list}

    df, name = cuffdiff_table_selector(fpkm_table, "epifactor_fpkm_Gianni", **chars)

    print df

    df.to_csv('./epifactor_Gianni.csv')

    x, y = df.shape

    factors = np.zeros(x)

    divider = (df.mean(axis=1) / 100).values

    for i in range(x):
        if divider[i] != 0:
            factors[i] = 1 / 1 #divider[i]
    # print factors


    import matplotlib.pyplot as plt

    for i in range(0, df.shape[0], 30):

        sub_df = df.ix[i:i+30, :]
        x, y = sub_df.shape

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_title("Notch Gene Expression")

        ## necessary variables
        ind = np.arange(x)  # the x locations for the groups
        width = 0.2  # the width of the bars

        rects1 = ax.bar(ind, sub_df['NEG_FPKM'].values * factors[i:i+30], width,
                        color='red')

        rects2 = ax.bar(ind + width, sub_df['LOW_FPKM'].values * factors[i:i+30], width,
                        color='green')

        rects3 = ax.bar(ind + 2 * width, sub_df['HIGH_FPKM'].values * factors[i:i+30], width,
                        color='blue')

        # rects4 = ax.bar(ind + 3 * width, df['HC'].values * factors, width,
        #                 color='yellow')

        # print rects1
        ax.set_xlim(-width, len(ind) + width)
        ax.set_ylabel('Expression')
        ax.set_title('Epigenetic factor in Fli-GFP part '+ str(i))
        xTickMarks = sub_df.index.values
        ax.set_xticks(ind + width)
        xtickNames = ax.set_xticklabels(xTickMarks)
        plt.setp(xtickNames, rotation=45, fontsize=5)

        ## add a legend
        ax.legend((rects1[0], rects2[0], rects3[0]), ('NEG', 'LOW', 'HIGH'), loc='best').draggable()
        plt.show()
        fig.savefig("./cook3/" + "histone"+str(i), dpi=200, facecolor='w', edgecolor='w',
                    orientation='portrait')
        plt.close('all')
    return df

def transform_scale(df):
    from sklearn.preprocessing import scale
    from pandas import DataFrame

    df = df.fillna(value=0)

    newdf = DataFrame(scale(df), index=df.index, columns=df.columns)

    return newdf


def bamTosam(path):
    file_names = os.listdir(path)

    for name in file_names:
        cmd = "samtools view -h -o " + name[:-4] + ".sam " + name
        os.system(cmd)


def HTSeq(path):
    file_names = [x for x in os.listdir(path) if x.endswith(".bam")]

    for name in file_names:
        cmd = "htseq-count -f bam -s no -m intersection-nonempty -i gene_id " + name + " /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.exon.anno.gtf > " + name[
                                                                                                                                                                      :-4] + "_count.txt"
        os.system(cmd)


def join_HTSeq_counts(path):
    file_names = [x for x in os.listdir(path) if x.endswith("_count.txt")]

    result_df = pd.read_csv(file_names[0], sep="\t", index_col=0, names=["gene_id", file_names[0][:file_names[0].find("_count.txt")]])

    for name in file_names[1:]:
        df = pd.read_csv(name, sep="\t", index_col=0, names=['gene_id', name[:name.find("_count.txt")]])

        result_df = pd.concat([result_df, df], axis=1, join='inner')

    result_df.to_csv("ec_vs_hsc_total.txt", sep="\t")

def common_elements(list1, list2):
    return list(set(list1) & set(list2))


def barplot_each_gene():
    f = open("./cook3/interesting_list", "r")
    #
    genes = [gene.strip().lower() for gene in f.readlines()]

    chars = {'gene_id': genes,
             'columns': ['NEG1_FPKM', 'NEG2_FPKM', 'LOW1_FPKM', 'LOW2_FPKM', 'HIGH1_FPKM', 'HIGH2_FPKM']}
    # #
    df = cuffdiff_table_selector("./cook3/genes.fpkm_tracking", '3-19', **chars)[0]

    # df.reset_index(level=0, inplace=True)
    new_df = pd.DataFrame(columns=['gene_id', 'FPKM', 'group', 'sd'])
    print df.columns.values
    for prefix in ['NEG', 'LOW', 'HIGH']:
        count = 0
        related_cols = []
        for column in df.columns.values:
            if column[:-6] == prefix:
                related_cols += [column]
                count += 1
        sub_df = df[related_cols].mean(axis=1).copy().to_frame()
        sub_df['group'] = [prefix for i in range(sub_df.shape[0])]
        sub_df['sd'] = df[related_cols].std(axis=1).values
        sub_df.reset_index(level=0, inplace=True)
        sub_df.columns = ['gene_id', 'FPKM', 'group', 'sd']
        new_df = pd.concat([new_df, sub_df])



if __name__ == "__main__":
    ###
    # chars = {'sample_1':'EC', 'sample_2':'HEC', 'q_value':0.05, 'log2(fold_change)': 0.0000000001}
    # cuffdiff_table_selector('./JEM/gene_exp.diff.tsv','EC_vs_HEC', **chars)
    #
    # chars = {'sample_1':'EC', 'sample_2':'HEC', 'q_value':0.05, 'log2(fold_change)': -0.0000000001}
    # cuffdiff_table_selector('./JEM/gene_exp.diff.tsv','EC_vs_HEC', **chars)
    # ###
    # chars = {'sample_1':'EC', 'sample_2':'HSC', 'q_value':0.05, 'log2(fold_change)': 0.0000000001}
    # cuffdiff_table_selector('./JEM/gene_exp.diff.tsv','EC_vs_HSC', **chars)
    #
    # chars = {'sample_1':'EC', 'sample_2':'HSC', 'q_value':0.05, 'log2(fold_change)': -0.0000000001}
    # cuffdiff_table_selector('./JEM/gene_exp.diff.tsv','EC_vs_HSC', **chars)
    # ###
    # chars = {'sample_1':'EC', 'sample_2':'HC', 'q_value':0.05, 'log2(fold_change)': 0.0000000001}
    # cuffdiff_table_selector('./JEM/gene_exp.diff.tsv','EC_vs_HC', **chars)
    #
    # chars = {'sample_1':'EC', 'sample_2':'HC', 'q_value':0.05, 'log2(fold_change)': -0.0000000001}
    # cuffdiff_table_selector('./JEM/gene_exp.diff.tsv','EC_vs_HC', **chars)

    # chars = {'sample_1':'HEC', 'sample_2':'HSC', 'q_value':0.05, 'log2(fold_change)': 0.0000000001}
    # cuffdiff_table_selector('./JEM/gene_exp.diff.tsv','HEC_vs_HSC', **chars)
    #
    # chars = {'sample_1':'HEC', 'sample_2':'HSC', 'q_value':0.05, 'log2(fold_change)': -0.0000000001}
    # cuffdiff_table_selector('./JEM/gene_exp.diff.tsv','HEC_vs_HSC', **chars)

    # f = open("./cook3/interesting_list", "r")
    # #
    # genes = [gene.strip().lower() for gene in f.readlines()]
    #
    # chars = {'gene_id': genes, 'columns':['NEG1_FPKM','NEG2_FPKM','LOW1_FPKM','LOW2_FPKM','HIGH1_FPKM','HIGH2_FPKM']}
    # # #
    # df = cuffdiff_table_selector("./cook3/genes.fpkm_tracking", '3-19', **chars)[0]
    #
    # # df.reset_index(level=0, inplace=True)
    # new_df = pd.DataFrame(columns=['gene_id', 'FPKM', 'group', 'sd'])
    # print df.columns.values
    # for prefix in ['NEG', 'LOW', 'HIGH']:
    #     count = 0
    #     related_cols = []
    #     for column in df.columns.values:
    #         if column[:-6] == prefix:
    #             related_cols += [column]
    #             count += 1
    #     sub_df = df[related_cols].mean(axis=1).copy().to_frame()
    #     sub_df['group'] = [prefix for i in range(sub_df.shape[0])]
    #     sub_df['sd'] = df[related_cols].std(axis=1).values
    #     sub_df.reset_index(level=0, inplace=True)
    #     sub_df.columns = ['gene_id', 'FPKM', 'group', 'sd']
    #     new_df = pd.concat([new_df, sub_df])
    #
    # new_df.to_csv('3-19.csv', index=False)
    #
    # chars = {'sample_1':'LOW', 'sample_2':'HIGH', 'q_value':0.05, 'log2(fold_change)': 1}
    # cuffdiff_table_selector('./cook3/gene_exp.diff','LOW_vs_HIGH', **chars)
    #
    # chars = {'sample_1':'LOW', 'sample_2':'HIGH', 'q_value':0.05, 'log2(fold_change)': -1}
    # cuffdiff_table_selector('./cook3/gene_exp.diff','LOW_vs_HIGH', **chars)

    # df = pd.read_csv("./JEM/JEM_total_quantile.txt", index_col=0, sep='\t')
    # print df
    #
    # df = df.ix[genes, ]
    # #
    # df.to_csv("JEM_ec_hec_q_0.05.tsv", sep='\t')


    # df = pd.read_csv('LOW_vs_HIGH_q_value_0.05_log2(fold_change)_1_sample_2_HIGH_sample_1_LOW.tsv', sep='\t')
    # df2 = pd.read_csv('LOW_vs_HIGH_up.tsv', sep='\t')
    #
    # genes1 = list(df['gene_id'].values)
    # genes2 = list(df2['gene_id'].values)
    # print len(genes1)
    # print len(genes2)
    # overlap = common_elements(genes1, genes2)
    #
    # f = open('down_in_HIGH_compare_NEGandLOW.txt', "w")
    # for g in overlap:
    #     f.write(g+"\n")
    # f.close()
    #


    # filenames = os.listdir("./cook3/david/")
    # path = "./cook3/david/"
    # #
    # for name in filenames:
    # #     os.system('mv '+path+name+' '+path+name[:-4]+"_david.tsv")
    #     try:
    #         df = pd.read_csv(path+name, sep='\t', index_col=0)
    #         del df['PValue']
    #         del df['List Total']
    #         del df['Pop Hits']
    #         del df['Pop Total']
    #         del df['Bonferroni']
    #         del df['FDR']
    #         df.to_csv(path+name, sep='\t')
    #     except:
    #         pass

    # df = pd.read_csv("./cook3/genes.fpkm_tracking", sep="\t", index_col=0)
    # genes = pd.read_csv("gene_list.txt")
    # genes = genes.ix[:, 0].values
    #
    # df = df.ix[genes, :]
    # df.to_csv("LOWHIGH.txt", sep="\t")

    # filenames = os.listdir("./cook3/david/")
    # path = "./cook3/david/"
    # #
    # for name in filenames:
    #     df = pd.read_csv(path + name, sep='\t', index_col=0)
    #     df = df.round({'Fold Enrichment':1})
    #     df.to_csv(path + name, sep='\t')

    # gene_exp_barplot('./cook3/genes.fpkm_tracking', './cook3/zebrafish_histone.txt')

    # notch_list_obj = open('./cook3/epifactorlist.txt', "r")
    #
    # notch_gene_list = [x.strip().lower() for x in notch_list_obj.readlines()]
    # notch_list_obj.close()
    #
    # chars = {"gene_id": notch_gene_list}
    #
    # df, name = cuffdiff_table_selector('./cook3/gene_exp.diff', "epifactor_Gianni", **chars)

    # print df

    df = pd.read_csv('./cook3/genes.fpkm_tracking', sep='\t')

    print df.mean(axis=0)
