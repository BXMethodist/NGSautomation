import os, pandas as pd

path = "./"

results = {}
results_left = {}
results_right={}
results_signle = {}
samples = [x for x in os.listdir(path) if os.path.isdir(x)]

if not path.endswith("/"):
    path += "/"
for sample in samples:
    file_path = path+sample+"/align_summary.txt"
    sample_obj = open(file_path, "r")
    info = sample_obj.readlines()
    sample_obj.close()

    sample_info = {}
    print sample
    for l in info:
        print l
        line = l.strip()
        if line.startswith("Left"):
            sample_info['left'] = {}
        elif line.startswith("Right"):
            sample_info['right'] = {}
        elif line.startswith("Aligned"):
            sample_info['all'] = {}
        elif line.startswith('Reads'):
            sample_info['read'] = {}

        if 'all' in sample_info:
            if line.startswith("Aligned pairs"):
                line_info = line.split(":")
                sample_info['all']['aligned'] = int(line_info[1])
            elif line.startswith("of these"):
                line_info = line.split(":")
                line_info = line_info[1].split()
                sample_info['all']['multiple'] = int(line_info[0])
            elif line.find("discordant") != -1:
                line_info = line.split()
                sample_info['all']['discordant'] = int(line_info[0])
            elif 'right' in sample_info:
                if line.startswith("Input"):
                    line_info = line.split(":")
                    sample_info['right']['input'] = int(line_info[1])
                elif line.startswith("Mapped"):
                    line_info = line.split(":")
                    line_info = line_info[1].split()
                    sample_info['right']['mapped'] = int(line_info[0])
                elif line.startswith("of these"):
                    line_info = line.split(":")
                    line_info = line_info[1].split()
                    sample_info['right']['multiple'] = int(line_info[0])
            elif 'left' in sample_info:
                if line.startswith("Input"):
                    line_info = line.split(":")
                    sample_info['left']['input'] = int(line_info[1])
                elif line.startswith("Mapped"):
                    line_info = line.split(":")
                    line_info = line_info[1].split()
                    line_info = line_info[1].split()
                    sample_info['read']['mapped'] = int(line_info[0])
                elif line.startswith("of these"):
                    line_info = line.split(":")
                    line_info = line_info[1].split()
                    sample_info['read']['multiple'] = int(line_info[0])

    if 'all' in sample_info:
        results[sample] = sample_info['all']
    if 'left' in sample_info:
        results_left[sample] = sample_info['left']
    if 'right' in sample_info:
        results_right[sample] = sample_info['right']
    if 'read' in sample_info:
        results_signle[sample] = sample_info['read']

df = pd.DataFrame(results)

df_left = pd.DataFrame(results_left)
df_right = pd.DataFrame(results_right)

df.to_csv('cook_pair_all.txt', sep='\t')
df_left.to_csv('cook_left.txt', sep='\t')
df_right.to_csv('cook_right.txt', sep='\t')



