import os

fast = ['SRR3472958',
'SRR3472951',
'SRR3472963',
'SRR3472948',
'SRR3472968']


path = '/home/tmhbxx3/archive/ec_vs_hsc/SRR/'

cmd = "/home/tmhbxx3/tools/sratoolkit/bin/fastq-dump -O ./FASTQ "

for f in fast:
    cur_cmd = cmd +path+f
    os.system(cur_cmd)