import os

error_files = [x for x in os.listdir('./') if x.find('.e') != -1]

errors = set()

failed = set()

for error_file in error_files:
    f = open('./' + error_file, "r")
    info = f.readlines()
    f.close()
    error = True

    for line in info:
        if line.startswith("# reads with at least one reported alignment:"):
            error = False
            percentage = float(line[line.find("(") + 1:line.find(")")-1])/100
        if line.startswith("# reads processed:"):
            reads_number = int(line.split()[-1])
        if line.startswith("# reads with alignments suppressed due to -m:"):
            ununique_reads = float(line[line.find("(") + 1:line.find(")")-1])/100
        if line.find("storage exhausted while writing file within file system module") != -1:
            errors.add(error_file[error_file.find('SRR') : error_file.find('.')])
            print error_file, "storage exhausted while writing file within file system module"
        elif line.find("fastq-dump.2.8.0") != -1:
            errors.add(error_file[error_file.find('SRR') : error_file.find('.')])
            print error_file, "fastq-dump.2.8.0"

    if error:
        failed.add(error_file[error_file.find('SRR') : error_file.find('.')])


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

print re_errors, len(re_errors), 'errors'
print re_failed, len(re_failed), 'failed'


[tmhbxx3@master report-5-30-2017]$ python bowtie_error.py
set(['SRR3591627', 'SRR1947844', 'SRR503400', 'SRR3090673', 'SRR1528849', 'SRR3178542', 'SRR611890', 'SRR3090674', 'SRR3090675', 'SRR3159917', 'SRR3126213', 'SRR3144849', 'SRR3144848', 'SRR3126214', 'SRR3591637', 'SRR3591625', 'SRR3144847', 'SRR3144846', 'SRR3144845', 'SRR3591617', 'SRR3591636', 'SRR3032859', 'SRR574819', 'SRR3032854', 'SRR3090670', 'SRR3090669', 'SRR3090668', 'SRR3363254', 'SRR3090667', 'SRR444466', 'SRR2079675', 'SRR3561326', 'SRR3561328', 'SRR3032848', 'SRR3032847', 'SRR3178537', 'SRR3032841', 'SRR3032842', 'SRR708094', 'SRR5060837', 'SRR488588', 'SRR3402852', 'SRR3471041', 'SRR3471042', 'SRR3098508', 'SRR3098504', 'SRR3098505', 'SRR3098503', 'SRR3032871', 'SRR3032872', 'SRR4437076', 'SRR4032239', 'SRR3363255', 'SRR3032877', 'SRR3032878', 'SRR3471048', 'SRR1947761', 'SRR3561331', 'SRR639066', 'SRR3624842', 'SRR3561332', 'SRR4420630', 'SRR1960157', 'SRR3098517', 'SRR3098516', 'SRR3098519', 'SRR3098518', 'SRR424619', 'SRR3032866', 'SRR3032865', 'SRR3561342', 'SRR3032860', 'SRR1515142', 'SRR1515140', 'SRR1515141', 'SRR4420629', 'SRR4420628', 'SRR424646', 'SRR524712', 'SRR3591616', 'SRR549368', 'SRR1960160', 'SRR3129109', 'SRR4032254', 'SRR3129106', 'SRR3129105', 'SRR3129104', 'SRR3471022', 'SRR488592', 'SRR3159928', 'SRR3159929', 'SRR3575341', 'SRR3471028', 'SRR3126158', 'SRR1515139', 'SRR3032834', 'SRR3090627', 'SRR3090628', 'SRR827530', 'SRR827532', 'SRR3575342', 'SRR1481512', 'SRR4032223', 'SRR3126161', 'SRR3126162', 'SRR3126163', 'SRR3126164', 'SRR708096', 'SRR3032853', 'SRR3471031', 'SRR3471033', 'SRR827521', 'SRR3471034', 'SRR3471037', 'SRR3471036', 'SRR3159930', 'SRR3032744', 'SRR3575339', 'SRR3180981', 'SRR5097107', 'SRR3180984', 'SRR3180985', 'SRR2155243', 'SRR2155242', 'SRR3129087', 'SRR3624832', 'SRR2155246', 'SRR2155245', 'SRR1599256', 'SRR4032230', 'SRR3561340', 'SRR3129082', 'SRR2079823', 'SRR1599254', 'SRR827520', 'SRR3178538', 'SRR3032758', 'SRR2969511', 'SRR3129081', 'SRR3032757', 'SRR1599253', 'SRR3591623', 'SRR5136487', 'SRR4299704', 'SRR4413986', 'SRR4299702', 'SRR4299703', 'SRR2155258', 'SRR3129095', 'SRR4032202', 'SRR4032203', 'SRR444442', 'SRR444443', 'SRR3591624', 'SRR2969521', 'SRR2969522', 'SRR3144854', 'SRR1599257', 'SRR3159915', 'SRR3591626', 'SRR3144850', 'SRR3144851', 'SRR3144852', 'SRR3144853', 'SRR3159918', 'SRR3591629', 'SRR3591628', 'SRR1599259']) 168