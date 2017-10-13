import os, pandas as pd

error_files = [x for x in os.listdir('./') if x.find('.e') != -1]

errors = set()

failed = set()

results = {}
for error_file in error_files:
    f = open('./' + error_file, "r")
    info = f.readlines()
    f.close()
    error = True

    name = error_file[error_file.find('SRR'): error_file.find('.e')]

    for line in info:
        if line.startswith("# reads with at least one reported alignment:"):
            error = False
            percentage = int(line[line.find(":") + 1:line.find("(")].strip())
        if line.startswith("# reads processed:"):
            reads_number = int(line.split()[-1])
        if line.startswith("# reads with alignments suppressed due to -m:"):
            ununique_reads = int(line[line.find(":") + 1:line.find("(")])
        if line.find("storage exhausted while writing file within file system module") != -1:
            errors.add(error_file[error_file.find('SRR') : error_file.find('.')])
        elif line.find("fastq-dump.2.8.0") != -1:
            errors.add(error_file[error_file.find('SRR') : error_file.find('.')])

    if error:
        failed.add(error_file[error_file.find('SRR') : error_file.find('.')])
    else:
        if name not in results:
            results[name] = [name, percentage, reads_number, ununique_reads]
        else:
            if int(results[name][1]) < int(percentage):
                results[name] = [name, percentage, reads_number, ununique_reads]

# print failed, len(failed)

print errors, len(errors)

re_errors = set()
re_failed = set()

together = failed.union(errors)

for fail in errors:
    duplicate = set()
    for error_file in error_files:
        if error_file.find(fail) != -1:
            duplicate.add(error_file)
    # print duplicate
    duplicate = sorted(list(duplicate), key=lambda x: int(x[x.find('.e')+2:]))
    # print duplicate

    f = open('./' + duplicate[-1], "r")
    info = f.readlines()
    f.close()
    error = True

    for line in info:
        if line.startswith("# reads with at least one reported alignment:"):
            error = False
            percentage = float(line[line.find("(") + 1:line.find(")") - 1]) / 100
        if line.startswith("# reads processed:"):
            reads_number = int(line.split()[-1])
        if line.startswith("# reads with alignments suppressed due to -m:"):
            ununique_reads = float(line[line.find("(") + 1:line.find(")") - 1]) / 100
        if line.find("storage exhausted while writing file within file system module") != -1:
            re_errors.add(fail)
        elif line.find("fastq-dump.2.8.0") != -1:
            re_errors.add(fail)

    if error:
        re_failed.add(error_file[error_file.find('SRR'): error_file.find('.')])

e_out = open('srr_errors.txt', 'w')
for e in re_errors:
    e_out.write(e+'\n')
e_out.close()

print re_errors, len(re_errors)
print re_failed, len(re_failed)

final_results = []

for value in results.values():
    final_results.append(value)

df = pd.DataFrame(final_results)
df.columns = ['sample_id', 'unique_mapped_reads', 'total_reads', 'ununique_mapped_reads']

df.to_csv('bowtie_srr_results.csv', index=None)

bowties = [b for b in os.listdir('../bowtie') if b.find('_') != -1 and b.endswith('.bowtie')]

underscore_bowtie = []

for b in bowties:
    name = b[:b.find('.bowtie')]
    names = name.split('_')

    total_reads = 0
    total_unique_mapped_reads = 0
    total_ununique_mapped_reads = 0

    for n in names:
        total_reads += results[n][2]
        total_unique_mapped_reads += results[n][1]
        total_ununique_mapped_reads += results[n][3]

    underscore_bowtie.append((name, total_unique_mapped_reads, total_reads, total_ununique_mapped_reads))

df = pd.DataFrame(underscore_bowtie)
df.columns = ['sample_id', 'unique_mapped_reads', 'total_reads', 'ununique_mapped_reads']

df.to_csv('bowtie_srr_with_combined_results.csv', index=None)



