import pandas as pd
import numpy as np
import urllib, re, sqlite3, csv


db_addr = "/home/tmhbxx3/archive/GEO_MetaDatabase/geoMetaData.db"
db = sqlite3.connect(db_addr)
db.text_factory = str

df = pd.read_csv("SREBP_samples.csv", sep=",", header=None)

gsmToSRR = []

groupToSRR = []

for i in range(df.shape[1]):
    group_SRRids = []
    gsmids = [x for x in df.ix[:, i].values if x is not np.nan]
    for gsmid in gsmids:
        query = db.execute('select SRA from GSM where GSM_ID = "'+gsmid+'"').fetchall()
        SRXlink = query[0][0]

        print SRXlink

        SRRids = []

        page = urllib.urlopen(SRXlink).read()

        SRRids += re.findall("SRR[0-9]+", page)
        SRRids += re.findall("DRR[0-9]+", page)
        SRRids += re.findall("ERX[0-9]+", page)

        group_SRRids += SRRids

        for SRRid in SRRids:
            gsmToSRR.append((gsmid, SRRid))
    groupToSRR.append(group_SRRids)

output = open("SREBP_single_cell_gsmToSRR.csv", "w")

for row in gsmToSRR:
    output.write(row[0]+"\t"+row[1]+"\n")

output.close()

output = open("SREBP_treat_SRR.csv", "w")
writer = csv.writer(output)
for row in groupToSRR:
    writer.writerow(row)

output.close()

db.close()