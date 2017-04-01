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