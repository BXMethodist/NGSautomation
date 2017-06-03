import os, pandas as pd

def HTSeq(path):
    file_names = [x for x in os.listdir(path) if x.endswith(".bam")]

    for name in file_names:
        cmd = "htseq-count -f bam -s no -m intersection-nonempty -i gene_id " + name + " /archive/tmhkxc48/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf > " + name[
                                                                                                                                                                      :-4] + "_count.txt"
        os.system(cmd)


def join_HTSeq_counts(path):
    file_names = [x for x in os.listdir(path) if x.endswith("_count.txt")]

    result_df = pd.read_csv(file_names[0], sep="\t", index_col=0, names=["gene_id", file_names[0][:file_names[0].find("_count.txt")]])

    for name in file_names[1:]:
        df = pd.read_csv(name, sep="\t", index_col=0, names=['gene_id', name[:name.find("_count.txt")]])

        result_df = pd.concat([result_df, df], axis=1, join='inner')

    result_df.to_csv("ec_vs_hsc_total.txt", sep="\t")

HTSeq('./')