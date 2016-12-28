import os

### tophat2 directory
### /home/tmhbxx3/tools/tophat-2.1.1/tophat
### /archive/tmhkxc48/ref_data/hg19/bowtie2/hg19

### fastq-dump
### /home/tmhbxx3/tools/sratoolkit/bin/fastq-dump

### hisat2
### /home/tmhbxx3/tools/hisat2-2.0.5/hisat2
### /home/tmhbxx3/archive/ref_data/hg19/hisat2/hg19


def runTopHat2(path, finishedjob=[], tophatIndex=" /archive/tmhkxc48/ref_data/hg19/bowtie2/hg19 "):
    cmd = "/home/tmhbxx3/tools/tophat-2.1.1/tophat2 --mate-std-dev 200 -p 8 -r 203 "

    if not path.endswith("/"):
        path += "/"

    tophatpath = path + "/tophat_output/"
    fastqpath = path+"FASTQ/"

    fileNames = os.listdir(fastqpath)

    if not os.path.isdir(tophatpath):
        os.system("mkdir "+tophatpath)

    for name in fileNames[200:]:
        if not name.endswith(".fastq"):
            continue

        if name[:name.find("_1.fastq")] in finishedjob:
            continue

        if name.endswith("_1.fastq"):
            reversePair = name[:name.find("_1.fastq")]+"_2.fastq"
            if reversePair in fileNames:
                outputPath = "-o " + tophatpath + name[:name.find("_1.fastq")]
                curCmd = cmd + outputPath + tophatIndex + fastqpath + name + " " + fastqpath + reversePair
                # print curCmd
                os.system(curCmd)
        elif name.endswith("_2.fastq"):
            pass
        else:
            outputPath = "-o ./tophat_output/" + name[:name.find(".fastq")]
            curCmd = cmd + outputPath + tophatIndex + fastqpath + name
            # print curCmd
            os.system(curCmd)
    return


def runFastqdump(path):
    if not os.path.isdir("./FASTQ"):
        os.system("mkdir FASTQ")

    cmd = "/home/tmhbxx3/tools/sratoolkit/bin/fastq-dump -O ./FASTQ --split-3 "

    if isinstance(path, list):
        for SRR in path:
            os.system(cmd + SRR)

    if isinstance(path, str):
        if not path.endswith("/"):
            path += "/"

        fileNames = os.listdir(path)
        for name in fileNames:
            if name.endswith(".sra"):
                os.system(cmd + path + name)
    return

def moveAndChangeNameForTophat(path):
    '''

    :param path: the directory contains tophat_output
    :return:
    '''

    if not path.endswith("/"):
        path += "/"

    input_path = path + "tophat_output/"
    output_path = path + "bamfiles/"

    if not os.path.isdir(output_path):
        os.system("mkdir "+output_path)

    SRRlist = os.listdir(input_path)

    for SRR_name in SRRlist:
        old_file_location = input_path+SRR_name+"/accepted_hits.bam"
        new_file_location = input_path+SRR_name+"/"+SRR_name+".bam"
        os.system("mv " + old_file_location + " "+ new_file_location)
        os.system("mv "+new_file_location+" "+output_path)
    return


def danposRecall(k, cutoff):
    wigPath = "/archive/tmhkxc48/BroadH3K4me3/broadpeak201401/H3K4me3/dregion/pooled/"

    wigFiles = [path for path in os.listdir(wigPath) if path.endswith("wig")]


    n = 0
    for wig in wigFiles:
        if 180*k <=n< 180*(k+1):
            cmd = "python /archive/tmhkxc48/tools/danposTemp/danpos.py dpeak "+wigPath+wig+" -q "+str(cutoff)+" -f 0 -z 0 -o /home/tmhbxx3/archive/KFH3K4me3/"+str(cutoff)+"cutoff"
            os.system(cmd)
        n+=1





if __name__ == "__main__":
    runTopHat2("/home/tmhbxx3/archive/SREBP", finishedjob=os.listdir("/home/tmhbxx3/archive/SREBP/tophat_output/")
               , tophatIndex=" /archive/tmhkxc48/ref_data/mm9/bowtie2/mm9 ")

    # moveAndChangeNameForTophat("/home/tmhbxx3/archive/SREBP")