import os

files = [x for x in os.listdir('./') if x.endswith('.bam')]
for bam in files:
    os.system('/home/tmhbxx3/archive/tools/bedtools2/bin/genomeCoverageBed -ibam '+bam+' -split -bg -g /archive/tmhkxc48/ref_data/mm9/mm9.chrom.sizes.xls >'+ bam[:-4]+'.bg')
    os.system('python /archive/tmhkxc48/lib/bedGraphLib.py nor2total ' + bam[:-4] + '.nor.bg 10000000000')
    os.system('/scratch/tmhdxz9/software/bin/bedGraphToBigWig ' + bam[:-4] + '.nor.bg /archive/tmhkxc48/ref_data/mm9/mm9.chrom.sizes.xls ' + bam[:-4] + '.nor.bw')