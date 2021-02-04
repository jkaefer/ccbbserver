
#!/usr/bin/env python
import sys
import csv
import json
import collections

# aliases
OrderedDict = collections.OrderedDict

src = '/tmp/data.log'
dst = '/tmp/data.json'
header = ['chr','pos','ref','alt','disease prob','FDR','transcript','ON/OFF splicing site','upstream intron length','exon length','downstream intron length','proximity acceptor site','proximity donor site','conservation','acceptor_site strength','donor_site strength','junction score change','PTM','Pfam']

data = []
with open(str(sys.argv[1]), 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='"')
    reader=list(reader)
    for row in range(len(reader)):
        if row==0:
            continue
        #if row[0].strip()[0] == '#':  #
        #    continue
        reader[row] = filter(None, reader[row])
        data.append(dict(OrderedDict(zip(header, reader[row]))))
      

print(data[0])
with open(str(sys.argv[2]), 'w') as jsonfile:
    json.dump(dict({"data":data}), jsonfile)

    

