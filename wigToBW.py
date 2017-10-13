import os

files = [x for x in os.listdir('./') if x.endswith('.wig')]

for wig in files:
    os.system('/home/tmhbxx3/archive/tools/ucsc/wigToBigWig -clip ' + wig + ' /archive/tmhkxc48/ref_data/hg19/hg19.chrom.sizes.xls '+wig[:-4]+'.bw')
    # os.system('/home/tmhbxx3/archive/tools/ucsc/wigToBigWig -clip ' + wig + ' /archive/tmhkxc48/ref_data/mm9/mm9.chrom.sizes.xls '+wig[:-4]+'.bw')
