/home/tmhbxx3/tools/tophat-2.1.1/tophat2 --mate-std-dev 200 -p 8 -r 203 -o /home/tmhbxx3/archive/SREBP/tophat_output/SRR2508836/  /archive/tmhkxc48/ref_data/mm9/bowtie2/mm9 /home/tmhbxx3/archive/SREBP/FASTQ/SRR2508836_1.fastq /home/tmhbxx3/archive/SREBP/FASTQ/SRR2508836_2.fastq



/home/tmhbxx3/tools/hisat2-2.0.5/hisat2-build /home/tmhbxx3/archive/ref_data/danRer10/danRer10.fa danRer10



## to count the reads length in fastq file
awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' CI5701_GCCAAT_L005_R1_001.fastq


track type=bigWig name="srebp2_bgsub" description="" visibility=2  db="mm9" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:20 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/srebp2_treat_GSM694314.bgsub.bw

track type=bigWig name="srebp2_sample" description="" visibility=2  db="mm9" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:20 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/GSM694314_srebp2.bw

track type=bigWig name="srebp2_input" description="" visibility=2  db="mm9" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:20 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/GSM694315_input.bw

track type=bigBed name="srebp2_input" description="" visibility=2  db="mm9" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:20 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/SREBP2/SREBP2_peaks.bb

cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/encode_ec_vs_hsc/cuffdiff_together -library-norm-method classic-fpkm --labels CD34+,HAoEC /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.exon.anno.gtf /home/tmhbxx3/archive/encode_ec_vs_hsc/bamfiles/SRR534325.bam /home/tmhbxx3/archive/encode_ec_vs_hsc/bamfiles/SRR545717.bam,/home/tmhbxx3/archive/encode_ec_vs_hsc/bamfiles/SRR545718.bam



tophat2 -p 8 -o ./tophat_out/CI5701 ~/archive/ref_data/danRer10/bowtie2/danRer10 CI5701_GCCAAT_L005_R2_001.fastq.gz
tophat2 -p 8 -o ./tophat_out/CI5703 ~/archive/ref_data/danRer10/bowtie2/danRer10 CI5703_CAGATC_L005_R2_001.fastq.gz
tophat2 -p 8 -o ./tophat_out/CI5705 ~/archive/ref_data/danRer10/bowtie2/danRer10 CI5705_ACTTGA_L005_R2_001.fastq.gz
tophat2 -p 8 -o ./tophat_out/CI5757 ~/archive/ref_data/danRer10/bowtie2/danRer10 CI5757_TAGCTT_L005_R2_001.fastq.gz
tophat2 -p 8 -o ./tophat_out/CI5758 ~/archive/ref_data/danRer10/bowtie2/danRer10 CI5758_GGCTAC_L005_R2_001.fastq.gz
tophat2 -p 8 -o ./tophat_out/CI5758 ~/archive/ref_data/danRer10/bowtie2/danRer10 CI5759_CTTGTA_L005_R2_001.fastq.gz


zfin  dre


/home/tmhbxx3/tools/bedtools2/bin/genomeCoverageBed -ibam NEG1.bam -split -bg -g /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt >CI5701.bg;
/home/tmhbxx3/tools/bedtools2/bin/genomeCoverageBed -ibam LOW1.bam -split -bg -g /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt >CI5703.bg;
/home/tmhbxx3/tools/bedtools2/bin/genomeCoverageBed -ibam HIGH1.bam -split -bg -g /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt >CI5705.bg;
/home/tmhbxx3/tools/bedtools2/bin/genomeCoverageBed -ibam NEG2.bam -split -bg -g /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt >CI5757.bg;
/home/tmhbxx3/tools/bedtools2/bin/genomeCoverageBed -ibam LOW2.bam -split -bg -g /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt >CI5758.bg;
/home/tmhbxx3/tools/bedtools2/bin/genomeCoverageBed -ibam HIGH2.bam -split -bg -g /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt >CI5759.bg;

python /archive/tmhkxc48/lib/bedGraphLib.py nor2total CI5701.bg 10000000000;
python /archive/tmhkxc48/lib/bedGraphLib.py nor2total CI5703.bg 10000000000;
python /archive/tmhkxc48/lib/bedGraphLib.py nor2total CI5705.bg 10000000000;
python /archive/tmhkxc48/lib/bedGraphLib.py nor2total CI5757.bg 10000000000;
python /archive/tmhkxc48/lib/bedGraphLib.py nor2total CI5758.bg 10000000000;
python /archive/tmhkxc48/lib/bedGraphLib.py nor2total CI5759.bg 10000000000;


/scratch/tmhdxz9/software/bin/bedGraphToBigWig CI5701.bg /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt CI5701.nor.bw;
/scratch/tmhdxz9/software/bin/bedGraphToBigWig CI5703.nor.bg /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt CI5703.nor.bw;
/scratch/tmhdxz9/software/bin/bedGraphToBigWig CI5705.nor.bg /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt CI5705.nor.bw;
/scratch/tmhdxz9/software/bin/bedGraphToBigWig CI5757.nor.bg /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt CI5757.nor.bw;
/scratch/tmhdxz9/software/bin/bedGraphToBigWig CI5758.nor.bg /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt CI5758.nor.bw;
/scratch/tmhdxz9/software/bin/bedGraphToBigWig CI5759.nor.bg /home/tmhbxx3/archive/ref_data/danRer10/chrominfo.txt CI5759.nor.bw;




track type=bigWig name="NEG1" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:1000 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5701.nor.bw


track type=bigWig name="NEG2" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:1000 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5757.nor.bw


track type=bigWig name="LOW1" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:1000 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5703.nor.bw


track type=bigWig name="LOW2" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:1000 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5758.nor.bw


track type=bigWig name="HIGH1" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:1000 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5705.nor.bw


track type=bigWig name="HIGH2" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:1000 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5759.nor.bw






#!/bin/bash
#PBS -r n
#PBS -N wig_spliter
#PBS -q mediummem
#PBS -o output.out
#PBS -e error.err
#PBS -m e
#PBS -M bxia@houstonmethodist.org
#PBS -l walltime=96:00:00
#PBS -l nodes=1:ppn=1
#PBS -l pmem=16000mb
cd /home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/profile/up

module load python/2.7.11
module load R/3.2.1

python /archive/tmhkxc48/tools/danposTemp/danpos.py profile /home/tmhbxx3/archive/SREBP/chipseq/wig/smooth_wig/srebp2_treat_GSM694314.bgsub.smooth.wig --genefile_paths /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.xls,/home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/EC_HEC_up_q0.05.txt_profilelist.txt,/home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/EC_HEC_down_q0.05.txt_profilelist.txt --wigfile_aliases SREBP2 --genefile_aliases all_genes,EC_vs_HEC_up_in_HEC,EC_vs_HEC_down_in_HEC --name /home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/EC_HEC --genomic_sites TSS,TTS,CSS,CTS,ESS,ETS --heatmap 1 --flank_up 3000 --flank_dn 3000

python /archive/tmhkxc48/tools/danposTemp/danpos.py profile /home/tmhbxx3/archive/SREBP/chipseq/wig/smooth_wig/srebp2_treat_GSM694314.bgsub.smooth.wig --genefile_paths /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.xls,/home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/EC_HEC_up_top500.txt_profilelist.txt,/home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/EC_HEC_down_top500.txt_profilelist.txt --wigfile_aliases SREBP2 --genefile_aliases all_genes,EC_vs_HEC_up_in_HEC_top500,EC_vs_HEC_down_in_HEC_top500 --name /home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/EC_HEC_top500 --genomic_sites TSS,TTS,CSS,CTS,ESS,ETS --heatmap 1 --flank_up 3000 --flank_dn 3000




python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dregion srebp2_treat_GSM694314.wig -b srebp2_control_GSM694315.wig --smooth_width 0 --frsz 200 --extend 200 --extend_dis 3000 --pheight 1e-8 -ep 1e-5 -o SREBP2_peaks




sort -k1,1 -k2,2n ./srebp2_treat_GSM694314.bgsub.Fnor.regions.xls|sed 1d |cut -f 1-3 >  ./sorted_regions.bed
~/tools/ucsc/bedToBigBed ./sorted_regions.bed /archive/tmhkxc48/ref_data/mm9/mm9.chrom.sizes.xls ./SREBP2_regions.bb



Hi Bo,
     This is your ID on cigwiki.houstonmethodist.org server.
     username:   boxia
     password:   123456
     Please change your password immediately.

You could upload track file to /var/www/html/trackhub/boxia .
The track link is http://cigwiki.houstonmethodist.org/trackhub/boxia/……
If you have any question, please tell me.


/home/tmhbxx3/archive/SREBP/chipseq/wig/SREBP2_peakswithc/pooled/srebp2_treat_GSM694314.bgsub.Fnor.wig

python /archive/tmhkxc48/tools/danposTemp/danpos.py selector /home/tmhbxx3/archive/SREBP/chipseq/wig/SREBP2_peakswithc/pooled/srebp2_treat_GSM694314.bgsub.Fnor.regions.xls --genicSelector TSS:-3000:3000 --gene_file notch.txt_profilelist.txt --gene_out ./selector/notch.xlsx --out ./selector/notch_pathway.xlsx

python /archive/tmhkxc48/tools/danposTemp/danpos.py selector /home/tmhbxx3/archive/SREBP/chipseq/wig/SREBP2_peakswithc/pooled/srebp2_treat_GSM694314.bgsub.Fnor.regions.xls --genicSelector TSS:-3000:3000 --gene_file EC_HEC_up.txt_profilelist.txt --gene_out ./selector/EC_HEC_up_gene.xlsx --out ./selector/EC_HEC_up_region.xlsx

python /archive/tmhkxc48/tools/danposTemp/danpos.py selector /home/tmhbxx3/archive/SREBP/chipseq/wig/SREBP2_peakswithc/pooled/srebp2_treat_GSM694314.bgsub.Fnor.regions.xls --genicSelector TSS:-3000:3000 --gene_file EC_HEC_down.txt_profilelist.txt --gene_out ./selector/EC_HEC_down_gene.xlsx --out ./selector/EC_HEC_down_region.xlsx

python /archive/tmhkxc48/tools/danposTemp/danpos.py selector /home/tmhbxx3/archive/SREBP/chipseq/wig/SREBP2_peakswithc/pooled/srebp2_treat_GSM694314.bgsub.Fnor.regions.xls --genicSelector TSS:-3000:3000 --gene_file notch_enriched.txt_profilelist.txt --gene_out ./selector/notch_enriched_gene.xlsx --out ./selector/notch_enriched_region.xlsx


python /archive/tmhkxc48/tools/danposTemp/danpos.py selector /home/tmhbxx3/archive/SREBP/chipseq/wig/SREBP2_peakswithc/pooled/srebp2_treat_GSM694314.bgsub.Fnor.regions.xls --genicSelector TSS:-3000:3000 --gene_file /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.xls --gene_out ./selector/all_gene.xlsx --out ./selector/all_region.xlsx

9
python /archive/tmhkxc48/tools/danposTemp/danpos.py profile /home/tmhbxx3/archive/SREBP/chipseq/wig/SREBP2_peakswithc/pooled/srebp2_treat_GSM694314.bgsub.Fnor.wig --genefile_paths /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.xls,/home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/EC_HEC_up.txt_profilelist.txt,/home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/EC_HEC_down.txt_profilelist.txt --wigfile_aliases SREBP2 --genefile_aliases all_genes,EC_HEC_up_in_HEC,EC_vs_HEC_down_in_HEC --genomic_sites TSS,TTS,CSS,CTS,ESS,ETS --heatmap 1 --flank_up 3000 --flank_dn 3000



python /archive/tmhkxc48/tools/danposTemp/danpos.py retrieveDNA /home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/selector/notch_enriched_region.xlsx /archive/tmhkxc48/ref_data/mm9/bowtie2/mm9.fa /home/tmhbxx3/archive/JEM_ec_hsc/chip-seq/notch_enriched_sequence.txt


fastx_trimmer -f 6 -i CI5701_GCCAAT_L005_R1_001.fastq -o CI5701_left_first_5_trimmed.fastq
fastx_trimmer -f 6 -i CI5703_CAGATC_L005_R1_001.fastq -o CI5703_left_first_5_trimmed.fastq
fastx_trimmer -f 6 -i CI5705_ACTTGA_L005_R1_001.fastq -o CI5705_left_first_5_trimmed.fastq
fastx_trimmer -f 6 -i CI5757_TAGCTT_L005_R1_001.fastq -o CI5757_left_first_5_trimmed.fastq
fastx_trimmer -f 6 -i CI5758_GGCTAC_L005_R1_001.fastq -o CI5758_left_first_5_trimmed.fastq
fastx_trimmer -f 6 -i CI5759_CTTGTA_L005_R1_001.fastq -o CI5759_left_first_5_trimmed.fastq

fastx_trimmer -f 6 -i CI5701_GCCAAT_L005_R2_001.fastq -o CI5701_right_first_5_trimmed.fastq
fastx_trimmer -f 6 -i CI5703_CAGATC_L005_R2_001.fastq -o CI5703_right_first_5_trimmed.fastq
fastx_trimmer -f 6 -i CI5705_ACTTGA_L005_R2_001.fastq -o CI5705_right_first_5_trimmed.fastq
fastx_trimmer -f 6 -i CI5757_TAGCTT_L005_R2_001.fastq -o CI5757_right_first_5_trimmed.fastq
fastx_trimmer -f 6 -i CI5758_GGCTAC_L005_R2_001.fastq -o CI5758_right_first_5_trimmed.fastq
fastx_trimmer -f 6 -i CI5759_CTTGTA_L005_R2_001.fastq -o CI5759_right_first_5_trimmed.fastq



#!/bin/bash
#PBS -r n
#PBS -N tophat
#PBS -q mediummem
#PBS -o output.out
#PBS -e error.err
#PBS -m e
#PBS -M bxia@houstonmethodist.org
#PBS -l walltime=96:00:00
#PBS -l nodes=1:ppn=8
#PBS -l pmem=16000mb

cd /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed

#tophat2 -p 8 -o ./tophat_out/NEG1 ~/archive/ref_data/danRer10/bowtie2/danRer10 /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5701_left_first_5_trimmed.fastq /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5701_right_first_5_trimmed.fastq

#tophat2 -p 8 -o ./tophat_out/LOW1 ~/archive/ref_data/danRer10/bowtie2/danRer10 /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5703_left_first_5_trimmed.fastq /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5703_right_first_5_trimmed.fastq

#tophat2 -p 8 -o ./tophat_out/HIGH1 ~/archive/ref_data/danRer10/bowtie2/danRer10 /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5705_left_first_5_trimmed.fastq /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5705_right_first_5_trimmed.fastq

#tophat2 -p 8 -o ./tophat_out/NEG2 ~/archive/ref_data/danRer10/bowtie2/danRer10 /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5757_left_first_5_trimmed.fastq /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5757_right_first_5_trimmed.fastq

#tophat2 -p 8 -o ./tophat_out/LOW2 ~/archive/ref_data/danRer10/bowtie2/danRer10 /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5758_left_first_5_trimmed.fastq /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5758_right_first_5_trimmed.fastq

#tophat2 -p 8 -o ./tophat_out/HIGH2 ~/archive/ref_data/danRer10/bowtie2/danRer10 /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5759_left_first_5_trimmed.fastq /home/tmhbxx3/archive/cooklab/Gianni/RNA-seq/FASTQ/trimmed/CI5759_right_first_5_trimmed.fastq

/home/tmhbxx3/scratch/GCF_OL


track type=bigWig name="NEG1" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5701.nor.bw

track type=bigWig name="NEG2" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5757.nor.bw

track type=bigWig name="LOW1" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5703.nor.bw

track type=bigWig name="LOW2" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5758.nor.bw

track type=bigWig name="HIGH1" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5705.nor.bw

track type=bigWig name="HIGH2" description="" visibility=2  db="danRer10" color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=off bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/Gianni/CI5759.nor.bw


python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dregion controlAR:doxAR -b controlAR:controlinput,doxAR:doxinput -u 1 --smooth_width 0 -c 25000000 --frsz 200 --extend 200 --extend_dis 3000 --pheight 1e-8 -ep 1e-5 -o compareclonalcut/dregion_controlAR_doxAR



tophat --mate-std-dev 200 -p 8 -r 200 -o Sample_B_kd1 /archive/tmhkxc48/ref_data/hg19/bowtie2/hg19 Sample_B_kd1/B_kd1_GCCAAT_L003_R1_001.fastq.gz Sample_B_kd1/B_kd1_GCCAAT_L003_R2_001.fastq.gz

ENCFF340BKZ.bowtie -b ENCFF647DBV.bowtie
ENCFF304NDL.bowtie -b ENCFF829GAK.bowtie
ENCFF580NZZ.bowtie -b ENCFF349FFX.bowtie


/home/tmhbxx3/archive/tools/ucsc/wigToBigWig -clip ./SRR385667.bgsub.Fnor.wig /archive/tmhkxc48/ref_data/hg19/hg19.chrom.sizes.xls ./SRR385667.bgsub.Fnor.bw;


track type=bigWig name="ENCFF002EEY" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/ENCFF002EEY.bgsub.Fnor.bw
track type=bigWig name="SRR1002521" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR1002521.bgsub.Fnor.bw
track type=bigWig name="SRR1151534" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR1151534.Fnor.bw
track type=bigWig name="SRR849201" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR849201.Fnor.bw
track type=bigWig name="ENCFF002EEX" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/ENCFF002EEX.bgsub.Fnor.bw
track type=bigWig name="ENCFF304NDL" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/ENCFF304NDL.bgsub.Fnor.bw
track type=bigWig name="ENCFF422BBZ" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/ENCFF422BBZ.bgsub.Fnor.bw
track type=bigWig name="SRR385664" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR385664.bgsub.Fnor.bw
track type=bigWig name="SRR385666" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR385666.bgsub.Fnor.bw
track type=bigWig name="SRR385667" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR385667.bgsub.Fnor.bw
track type=bigWig name="SRR385669" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR385669.bgsub.Fnor.bw
track type=bigWig name="SRR501051" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR501051.bgsub.Fnor.bw
track type=bigWig name="SRR501055" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR501055.bgsub.Fnor.bw
track type=bigWig name="SRR501059" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR501059.bgsub.Fnor.bw
track type=bigWig name="SRR501062" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR501062.bgsub.Fnor.bw
track type=bigWig name="SRR849202" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR849202.Fnor.bw

track type=bigWig name="SRR1550922" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR1550922.bgsub.Fnor.bw
track type=bigWig name="SRR385660" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR385660.Fnor.bw
track type=bigWig name="SRR385661" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR385661.Fnor.bw
track type=bigWig name="SRR2005769" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR2005769.bgsub.Fnor.bw
track type=bigWig name="SRR2571498" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR2571498.bgsub.Fnor.bw
track type=bigWig name="SRR2571502" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR2571502.bgsub.Fnor.bw
track type=bigWig name="SRR2571506" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR2571506.bgsub.Fnor.bw
track type=bigWig name="SRR2571510" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR2571510.bgsub.Fnor.bw
track type=bigWig name="SRR2571514" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR2571514.bgsub.Fnor.bw
track type=bigWig name="SRR2571528" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR2571528.bgsub.Fnor.bw
track type=bigWig name="SRR385658" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR385658.bgsub.Fnor.bw
track type=bigWig name="SRR385659" description="" visibility=2 db=mm9 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/BMI1/SRR385659.bgsub.Fnor.bw


/54Pplus_S35_L006_/align_summary.txt
Left reads:

          Input     :  18359040

           Mapped   :    266088 ( 1.4% of input)

            of these:     67318 (25.3%) have multiple alignments (3409 have >20)

Right reads:

          Input     :  18359040

           Mapped   :    163771 ( 0.9% of input)

            of these:     30078 (18.4%) have multiple alignments (3384 have >20)

 1.2% overall read mapping rate.



Aligned pairs:    134345

     of these:     18596 (13.8%) have multiple alignments

                    6925 ( 5.2%) are discordant alignments

 0.7% concordant pair alignment rate.

./T4tdplusI_S29_L005_/align_summary.txt
Left reads:

          Input     :  40850695

           Mapped   :     26528 ( 0.1% of input)

            of these:      4360 (16.4%) have multiple alignments (646 have >20)

Right reads:

          Input     :  40850695

           Mapped   :     21114 ( 0.1% of input)

            of these:      2766 (13.1%) have multiple alignments (639 have >20)

 0.1% overall read mapping rate.



Aligned pairs:     19286

     of these:      2216 (11.5%) have multiple alignments

                    1135 ( 5.9%) are discordant alignments

 0.0% concordant pair alignment rate.

./T54tdplusI_S31_L005_/align_summary.txt
Left reads:

          Input     :  34089900

           Mapped   :     47998 ( 0.1% of input)

            of these:      4466 ( 9.3%) have multiple alignments (1120 have >20)

Right reads:

          Input     :  34089900

           Mapped   :     44529 ( 0.1% of input)

            of these:      3936 ( 8.8%) have multiple alignments (1119 have >20)

 0.1% overall read mapping rate.



Aligned pairs:     42062

     of these:      3460 ( 8.2%) have multiple alignments

                    2087 ( 5.0%) are discordant alignments

 0.1% concordant pair alignment rate.

./9Iplus_S37_L006_/align_summary.txt
Left reads:

          Input     :  26840742

           Mapped   :    163891 ( 0.6% of input)

            of these:     16850 (10.3%) have multiple alignments (2925 have >20)

Right reads:

          Input     :  26840742

           Mapped   :    136192 ( 0.5% of input)

            of these:     10230 ( 7.5%) have multiple alignments (2924 have >20)

 0.6% overall read mapping rate.



Aligned pairs:    127552

     of these:      8579 ( 6.7%) have multiple alignments

                    4956 ( 3.9%) are discordant alignments

 0.5% concordant pair alignment rate.

./T4td-P_S24_L005_/align_summary.txt
Left reads:

          Input     :  41207091

           Mapped   :    180583 ( 0.4% of input)

            of these:     48154 (26.7%) have multiple alignments (3946 have >20)

Right reads:

          Input     :  41207091

           Mapped   :    144885 ( 0.4% of input)

            of these:     34803 (24.0%) have multiple alignments (3868 have >20)

 0.4% overall read mapping rate.



Aligned pairs:    129348

     of these:     29096 (22.5%) have multiple alignments

                    5341 ( 4.1%) are discordant alignments

 0.3% concordant pair alignment rate.

./T54td-I_S30_L005_/align_summary.txt
Left reads:

          Input     :  31047262

           Mapped   :    162428 ( 0.5% of input)

            of these:     26722 (16.5%) have multiple alignments (2466 have >20)

Right reads:

          Input     :  31047262

           Mapped   :    131641 ( 0.4% of input)

            of these:     18881 (14.3%) have multiple alignments (2466 have >20)

 0.5% overall read mapping rate.



Aligned pairs:    117485

     of these:     15767 (13.4%) have multiple alignments

                    7530 ( 6.4%) are discordant alignments

 0.4% concordant pair alignment rate.

./4I-_S38_L006_/align_summary.txt
Left reads:

          Input     :  26590627

           Mapped   :   1680590 ( 6.3% of input)

            of these:    515039 (30.6%) have multiple alignments (9645 have >20)

Right reads:

          Input     :  26590627

           Mapped   :   1226884 ( 4.6% of input)

            of these:    359875 (29.3%) have multiple alignments (9645 have >20)

 5.5% overall read mapping rate.



Aligned pairs:   1023007

     of these:    299500 (29.3%) have multiple alignments

                   54255 ( 5.3%) are discordant alignments

 3.6% concordant pair alignment rate.

./4Iplus_S39_L006_/align_summary.txt
Left reads:

          Input     :  27811391

           Mapped   :   1536371 ( 5.5% of input)

            of these:    431117 (28.1%) have multiple alignments (6471 have >20)

Right reads:

          Input     :  27811391

           Mapped   :    988517 ( 3.6% of input)

            of these:    243602 (24.6%) have multiple alignments (6440 have >20)

 4.5% overall read mapping rate.



Aligned pairs:    793373

     of these:    178875 (22.5%) have multiple alignments

                   15956 ( 2.0%) are discordant alignments

 2.8% concordant pair alignment rate.

./4Pplus_S33_L006_/align_summary.txt
Left reads:

          Input     :  20921543

           Mapped   :   2238173 (10.7% of input)

            of these:    455308 (20.3%) have multiple alignments (57913 have >20)

Right reads:

          Input     :  20921543

           Mapped   :   1750564 ( 8.4% of input)

            of these:    311882 (17.8%) have multiple alignments (57642 have >20)

 9.5% overall read mapping rate.



Aligned pairs:   1592750

     of these:    250706 (15.7%) have multiple alignments

                   89838 ( 5.6%) are discordant alignments

 7.2% concordant pair alignment rate.

./T4tdplusP_S25_L005_/align_summary.txt
Left reads:

          Input     :  34515371

           Mapped   :   6167055 (17.9% of input)

            of these:   2039898 (33.1%) have multiple alignments (85031 have >20)

Right reads:

          Input     :  34515371

           Mapped   :   4669687 (13.5% of input)

            of these:   1448317 (31.0%) have multiple alignments (82184 have >20)

15.7% overall read mapping rate.



Aligned pairs:   3763478

     of these:   1102078 (29.3%) have multiple alignments

                   92292 ( 2.5%) are discordant alignments

10.6% concordant pair alignment rate.

./T54tdplusP_S27_L005_/align_summary.txt
Left reads:

          Input     :  32863341

           Mapped   :   1499498 ( 4.6% of input)

            of these:    337430 (22.5%) have multiple alignments (23927 have >20)

Right reads:

          Input     :  32863341

           Mapped   :   1220751 ( 3.7% of input)

            of these:    250208 (20.5%) have multiple alignments (23509 have >20)

 4.1% overall read mapping rate.



Aligned pairs:   1067783

     of these:    203627 (19.1%) have multiple alignments

                   45149 ( 4.2%) are discordant alignments

 3.1% concordant pair alignment rate.

./9I-_S36_L006_/align_summary.txt
Left reads:

          Input     :  31653276

           Mapped   :    132728 ( 0.4% of input)

            of these:     32427 (24.4%) have multiple alignments (1705 have >20)

Right reads:

          Input     :  31653276

           Mapped   :     81250 ( 0.3% of input)

            of these:     11675 (14.4%) have multiple alignments (1705 have >20)

 0.3% overall read mapping rate.



Aligned pairs:     70995

     of these:      7629 (10.7%) have multiple alignments

                    5387 ( 7.6%) are discordant alignments

 0.2% concordant pair alignment rate.

./T4td-I_S28_L005_/align_summary.txt
Left reads:

          Input     :  34835784

           Mapped   :   3033544 ( 8.7% of input)

            of these:    466176 (15.4%) have multiple alignments (60690 have >20)

Right reads:

          Input     :  34835784

           Mapped   :   2601181 ( 7.5% of input)

            of these:    359512 (13.8%) have multiple alignments (60599 have >20)

 8.1% overall read mapping rate.



Aligned pairs:   2370881

     of these:    300268 (12.7%) have multiple alignments

                  146860 ( 6.2%) are discordant alignments

 6.4% concordant pair alignment rate.

./T54td-P_S26_L005_/align_summary.txt
Left reads:

          Input     :  40176717

           Mapped   :   1159794 ( 2.9% of input)

            of these:    213978 (18.4%) have multiple alignments (3515 have >20)

Right reads:

          Input     :  40176717

           Mapped   :    896296 ( 2.2% of input)

            of these:    147036 (16.4%) have multiple alignments (3472 have >20)

 2.6% overall read mapping rate.



Aligned pairs:    761416

     of these:    116053 (15.2%) have multiple alignments

                   35049 ( 4.6%) are discordant alignments

 1.8% concordant pair alignment rate.

./4P-_S32_L006_/align_summary.txt
Left reads:

          Input     :  33028664

           Mapped   :   1207195 ( 3.7% of input)

            of these:    352171 (29.2%) have multiple alignments (24666 have >20)

Right reads:

          Input     :  33028664

           Mapped   :    915008 ( 2.8% of input)

            of these:    239620 (26.2%) have multiple alignments (24413 have >20)

 3.2% overall read mapping rate.



Aligned pairs:    805705

     of these:    189287 (23.5%) have multiple alignments

                   31469 ( 3.9%) are discordant alignments

 2.3% concordant pair alignment rate.

./54P-_S34_L006_/align_summary.txt
Left reads:

          Input     :  26709917

           Mapped   :   1539919 ( 5.8% of input)

            of these:    177287 (11.5%) have multiple alignments (17991 have >20)

Right reads:

          Input     :  26709917

           Mapped   :   1262379 ( 4.7% of input)

            of these:    114336 ( 9.1%) have multiple alignments (17976 have >20)

 5.2% overall read mapping rate.



Aligned pairs:   1169591

     of these:     84703 ( 7.2%) have multiple alignments

                   39080 ( 3.3%) are discordant alignments

 4.2% concordant pair alignment rate.


 tophat2 -p 8 --mate-std-dev 200 -r 200 -o ./T4tdplusP_S25_L005_ -G /archive/tmhkxc48/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf -T /archive/tmhkxc48/ref_data/mm9/bowtie2/mm9  T4tdplusP_S25_L005_R1_001.fastq T4tdplusP_S25_L005_R2_001.fastq


cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/TERT/TERT2/bam -library-norm-method classic-fpkm --labels TERTPlus,TERTMinus /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf TERTplus.bam TERTminus.bam


/home/tmhbxx3/archive/tools/ucsc/wigToBigWig -clip ./GSM947525.Fnor.wig /archive/tmhkxc48/ref_data/hg19/hg19.chrom.sizes.xls ./GSM947525.Fnor.bw;


track type=bigWig name="GSM1541011.bgsub.Fnor.bw" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/GSM1541011.bgsub.Fnor.bw


track type=bigWig name="GSM1541012.bgsub.Fnor.bw" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/GSM1541012.bgsub.Fnor.bw


track type=bigWig name="GSM947525.Fnor.bw" description="" visibility=2 db=hg19 color=0,0,255 maxHeightPixels=30:30:30 windowingFunction=maximum viewLimits=0:500 autoScale=auto bigDataUrl=http://cigwiki.houstonmethodist.org/trackhub/boxia/GSM947525.Fnor.bw


cuffdiff -p 8 --output-dir /home/tmhbxx3/archive/JQ1/GSE84520_macrophage/RNA-seq/cuffdiff -library-norm-method classic-fpkm --labels M0_UT,M1_IFN,M2_IL4 /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf SRR3929129.bam,SRR3929136.bam,SRR3929143.bam SRR3929134.bam,SRR3929141.bam,SRR3929148.bam SRR3929133.bam,SRR3929140.bam,SRR3929147.bam


cuffdiff -p 8 --output-dir /home/tmhbxx3/scratch/CIG_lg/RNA-seq/cuffdiff_MSC -library-norm-method classic-fpkm --labels MSC,MSC /archive/tmhkxc48/ref_data/hg19/hg19.ucscgenes.knowngene.exon.anno.gtf MSC.bam MSC.bam



