import os

path = "/home/tmhbxx3/archive/SREBP"

if not path.endswith("/"):
    path += "/"

output_dir = path + "cuffdiff_norm_classic-fpkm"

cmd = "cuffdiff -p 9 --output-dir "+output_dir

cmd += " -library-norm-method classic-fpkm --labels EC,FL,T1,T2 /home/tmhbxx3/archive/ref_data/mm9/mm9.20150218.knownGene.exon.anno.gtf "

if not os.path.isdir(output_dir):
    os.system("mkdir " + output_dir)

input_file = open(path+"/code/SREBP_treat_SRR.csv", "r")

for line in input_file.readlines():
    line = line.rstrip()
    cmd += line + " "

os.system(cmd)