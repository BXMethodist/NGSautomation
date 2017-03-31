#This module is used to convert gene id between different sources, transcript ID and gene_ID

import numpy as np, pandas as pd, os

def transcriptID_to_gene_ID(source, from_type, to_type, from_id, to_id, gtf):
    if source == "UCSC":
        return UCSC_transID_geneID(from_type, to_type, from_id, to_id, gtf)
    if source == 'Ensemble':
        return Ensemble_transID_geneID(from_type, to_type, from_id, to_id, gtf)
    if source == 'RefSeq':
        return RefSeq_transID_geneID(from_type, to_type, from_id, to_id, gtf)

def UCSC_transID_geneID(from_type, to_type, from_id, to_id, gtf):
    """
    :param from_type: transcript or gene
    :param to_type: transcript or gene
    :param from_id: a file containing a list of ids
    :param to_id: a file containing a list of ids
    :return:
    """
    # f = open(from_id, 'r')
    # info = f.readlines()
    # f.close()
    #
    # from_id = set()
    # for line in info:
    #     line = line.strip()
    #     from_id.add(line)

    gtf = open(gtf, 'r')
    gtf_info = gtf.readlines()
    gtf.close()

    gene_transcript_map = {}
    transcript_gene_map = {}

    for line in gtf_info:
        line = line.split('\t')
        gene, transcript = line[-1].split(";")[0:2]
        # print gene, transcript
        # break
        gene_transcript_map[gene[gene.find('\"')+1:gene.rfind('\"')]] = \
            transcript[transcript.find('\"')+1:transcript.rfind('\"')]
        transcript_gene_map[transcript[transcript.find('\"')+1:transcript.rfind('\"')]] = \
            gene[gene.find('\"') + 1:gene.rfind('\"')]

    result_ids = set()
    # print transcript_gene_map
    if from_type == 'transcript':
        for id in from_id:
            if id in transcript_gene_map:
                result_ids.add(transcript_gene_map[id])
            else:
                print id
    else:
        for id in from_id:
            if id in gene_transcript_map:
                result_ids.add(gene_transcript_map[id])
            else:
                print id

    result_file = open('UCSC'+to_type+'.txt', "w")
    for id in result_ids:
        result_file.write(id + '\n')
    result_file.close()
    return result_ids

def Ensemble_transID_geneID(from_type, to_type, from_id, to_id, gtf):
    """
    :param from_type: transcript or gene
    :param to_type: transcript or gene
    :param from_id: a file containing a list of ids
    :param to_id: a file containing a list of ids
    :return:
    """
    pass

def RefSeq_transID_geneID(from_type, to_type, from_id, to_id, gtf):
    """
    :param from_type: transcript or gene
    :param to_type: transcript or gene
    :param from_id: a file containing a list of ids
    :param to_id: a file containing a list of ids
    :return:
    """
    pass


df = pd.read_csv("SREBP_binding_gene.txt", sep='\t', index_col=0)

from_id = df.index.tolist()


transcriptID_to_gene_ID('UCSC',
                        'transcript',
                        'gene',
                        from_id, '',
                        '/archive/tmhkxc48/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf')

