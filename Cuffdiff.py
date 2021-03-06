from itertools import combinations
import os, pandas as pd

combs = [(x, y) for x in combinations([1,2,3], 2) for y in combinations([4,5,6],2)]

# for x, y in combs:
#     f = open('_'.join([str(i) for i in x]) + '_'.join([str(j) for j in y])+'.pbs', 'w')
#     f.write('#!/bin/bash\n')
#     f.write('#PBS -r n\n')
#     f.write('#PBS -m e\n')
#     f.write('#PBS -M bxia@houstonmethodist.org\n')
#     f.write('#PBS -l walltime=96:00:00\n')
#     f.write('#PBS -l nodes=1:ppn=8\n')
#     f.write('#PBS -l pmem=16000mb\n')
#     f.write('#PBS -q mediummem\n')
#     f.write('module load python/2.7.11\n')
#     f.write('cd /home/tmhbxx3/archive/PAF/bams\n')
#     f.write('cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/PAF/bams/'+ '_'.join([str(i) for i in x]) + '_'.join([str(j) for j in y]) + ' -library-norm-method classic-fpkm --labels WT,PAF /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf '+ ','.join([str(i) + '.bam' for i in x]) + ' ' +  ','.join([str(j)+'.bam' for j in y]) + '\n')
#     f.close()
#     os.system('qsub ' + '_'.join([str(i) for i in x]) + '_'.join([str(j) for j in y])+'.pbs')


for comb in combs:
    x, y = comb
    folder = './' + '_'.join([str(i) for i in x]) + '_'.join([str(j) for j in y]) +'/'
    df = pd.read_csv(folder+'gene_exp.diff', sep='\t')
    df = df[df['significant'] == 'yes']
    print folder, df.shape



path = "/home/tmhbxx3/archive/encode_ec_vs_hsc/bamfiles/"

# if not path.endswith("/"):
#     path += "/"
#
# bamfiles = path+"bamfiles/"
#
# file_names = [x for x in os.listdir(bamfiles) if x.endswith(".bam")]
#
# for name in file_names:
#     output_dir = path + "cuffdiff/" + name[:-4]
#
#     if not os.path.isdir(output_dir):
#         os.system("mkdir " + output_dir)
#
#     cmd = "cuffdiff -p 8 --output-dir "+output_dir
#
#     cmd += " -library-norm-method classic-fpkm  --labels "+ name[:-4] +","+name[:-4]+ " /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.exon.anno.gtf "+bamfiles+name+" "+bamfiles+name
#
#     os.system(cmd)

cmd = "cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/encode_ec_vs_hsc/cuffdiff_together -library-norm-method classic-fpkm " \
      "--labels CD34+,HAoEC /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.exon.anno.gtf SRR534325.bam SRR545717.bam,SRR545718.bam"



# cmd += " -library-norm-method classic-fpkm --labels NEG,Low,High /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf "
# --labels NEG1,NEG2,Low1,Low2,High1,High2

# cmd += " -library-norm-method classic-fpkm --labels NEG1,NEG2,Low1,Low2,High1,High2 /home/tmhbxx3/archive/ref_data/danRer10/danRer10.ensembl_common_name.gtf "


# input_file = ["CI5701_GCCAAT_L005_", "CI5757_TAGCTT_L005_",
#               "CI5703_CAGATC_L005_", "CI5758_GGCTAC_L005_",
#               "CI5705_ACTTGA_L005_", "CI5759_CTTGTA_L005_"]
#
#
# for line in input_file:
#     line = line.rstrip()
#     samples = line.split(",")
#     for sample in samples:
#         cmd += bamfiles+sample+".bam,"
#     cmd = cmd[:-1]
#     cmd += " "
#
# print cmd
# os.system(cmd)

import os

path = "/home/tmhbxx3/archive/JEM_ec_hsc/bamfiles"

names = [x for x in os.listdir(path) if x.endswith(".bam")]

cmd = "cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/JEM_ec_hsc/cuffdiff_together -library-norm-method classic-fpkm " \
      "--labels EC,HEC,HSC,HC /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf SRR1653121.bam," \
      "SRR1653122.bam,SRR1653123.bam,SRR1653124.bam SRR1653125.bam,SRR1653126.bam,SRR1653127.bam,SRR1653128.bam " \
      "SRR1653129.bam,SRR1653130.bam,SRR1653131.bam,SRR1653132.bam SRR1653133.bam,SRR1653134.bam,SRR1653135.bam,SRR1653136.bam"

os.system(cmd)

cmd = "cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/JEM_ec_hsc/cuffdiff_individual -library-norm-method classic-fpkm " \
      "--labels EC1,EC2,EC3,EC4,HEC1,HEC2,HEC3,HEC4,HSC1,HSC2,HSC3,HSC4,HC1,HC2,HC3,HC4 /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf SRR1653121.bam SRR1653122.bam SRR1653123.bam SRR1653124.bam SRR1653125.bam SRR1653126.bam SRR1653127.bam SRR1653128.bam SRR1653129.bam SRR1653130.bam SRR1653131.bam SRR1653132.bam SRR1653133.bam SRR1653134.bam SRR1653135.bam SRR1653136.bam"

os.system(cmd)

import os

cmd = "cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/cuffdiff_individual -library-norm-method classic-fpkm --labels NEG1,LOW1,HIGH1,NEG2,LOW2,HIGH2 /home/tmhbxx3/archive/ref_data/danRer10/danRer10.ensembl_common_name.gtf NEG1.bam LOW1.bam HIGH1.bam NEG2.bam LOW2.bam HIGH2.bam"

os.system(cmd)


cmd = "cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/cuffdiff_together -library-norm-method classic-fpkm --labels NEG,LOW,HIGH /home/tmhbxx3/archive/ref_data/danRer10/danRer10.ensembl_common_name.gtf NEG1.bam,NEG2.bam LOW1.bam,LOW2.bam HIGH1.bam,HIGH2.bam"
os.system(cmd)


Erin

cmd = "cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/Erin/RNA-seq/cuffdiff/individual -library-norm-method classic-fpkm --labels CAA,CAB,CAC,CSA,CSB,CSC,WAA,WAB,WAC,WSA,WSB,WSC /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf CAA_.bam,CAB_.bam,CAC_.bam,CSA_.bam,CSB_.bam,CSC_.bam,WAA_.bam,WAB_.bam,WAC_.bam,WSA_.bam,WSB_.bam,WSC_.bam"

PAF

cmd = "cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/PAF/cdts-wh.genomics.cn/F17FTSUSAT0305_MUShkkN/IGV/bam/cuffdiff -library-norm-method classic-fpkm --labels WT,PAF /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf 1.sorted.bam,2.sorted.bam,3.sorted.bam 4.sorted.bam,5.sorted.bam,6.sorted.bam"