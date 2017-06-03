import os, pandas as pd

def wigToBigWig(wigfilePath, wigToBigWigPath='/home/tmhbxx3/archive/tools/ucsc/wigToBigWig', chromsize=''):
    outputname = wigfilePath[:-4]+".bw"
    print " ".join([wigToBigWigPath, '-clip', wigfilePath, chromsize, outputname])
    # os.system(" ".join([wigToBigWigPath, wigfilePath, chromsize, outputname]))

def wigToBigWig_main(search, metadata, wigfilesPath):
    """
    :param search: the result from chipseqpair search
    :param metadata: the result from ChipSeqpair query
    :param wigfilesPath: the folder path containing all the wig files
    :return:
    """
    search_df = pd.read_csv(search, index_col=0)
    metadata_df = pd.read_csv(metadata, index_col=0, sep='\t')

    wigfilesPath = wigfilesPath + '/' if not wigfilesPath.endswith('/') else wigfilesPath + ''

    wigs = [wigfilesPath+x for x in os.listdir(wigfilesPath) if x.endswith('.wig')]

    for wig in wigs:
        if wig.find('ENC') != -1:
            name = wig[wig.rfind('/') + 1:wig.rfind('/') + 1 + wig[wig.rfind('/') + 1:].find('.')]
            species = search_df.ix[name, 'Organism']
            if species == 'Homo sapiens':
                chromsize = '/archive/tmhkxc48/ref_data/hg19/hg19.chrom.sizes.xls'
            elif species == "Mus musculus":
                chromsize = '/archive/tmhkxc48/ref_data/mm9/mm9.chrom.sizes.xls'
            wigToBigWig(wig, chromsize=chromsize)
        elif wig.find('SRR')!= -1:
            name = wig[wig.rfind('/') + 1:wig.rfind('/') + 1+wig[wig.rfind('/') + 1:].find('.')]
            GSM_ID = metadata_df.ix[name, 'GSM_ID']
            species = search_df.ix[GSM_ID, 'Organism']

            if species == 'Homo sapiens':
                chromsize = '/archive/tmhkxc48/ref_data/hg19/hg19.chrom.sizes.xls'
            elif species == "Mus musculus":
                chromsize = '/archive/tmhkxc48/ref_data/mm9/mm9.chrom.sizes.xls'
            wigToBigWig(wig, chromsize=chromsize)
        else:
            print wig

# wigToBigWig_main('Search_ResultHomo_sapiensWithBMI1.csv', 'BMI1_chipseq_metadata.txt', './')

files = [x for x in os.listdir('.') if x.endswith('.wig')]

for wigfilePath in files:
    wigToBigWig(wigfilePath, wigToBigWigPath='/home/tmhbxx3/archive/tools/ucsc/wigToBigWig', chromsize='/archive/tmhkxc48/ref_data/mm9/mm9.chrom.sizes.xls')