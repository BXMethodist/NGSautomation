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
