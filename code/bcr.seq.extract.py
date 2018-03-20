import csv
import argparse
from math import log as ln



def p(n, N):
    """ Relative abundance """
    if n is  0:
        return 0
    else:
        return (float(n)/N) * ln(float(n)/N)

def sdi(data):
    N = sum(data.values())

    return -sum(p(n, N) for n in data.values() if n is not 0)

#============================
ap = argparse.ArgumentParser()
ap.add_argument('input', help='<fastq.gz> file')
ap.add_argument('out_prefix', help='Directory to be created to save <fastq> files. Each <fastq> file is obtained based on the barcode')
args = ap.parse_args()

#ESO_1_sorted_inf        82.338987       81220   CASSYVGNTGELFF  CTGCTGTCGGCTGCTCCCTCCCAGACATCTGTGTACTTCTGTGCCAGCAGTTACGTCGGGAACACCGGGGAGCTGTTTTTTGGAGAA TCRBV06-05*01   TCRBD02-01      TCRBJ02-02*01

set_cdr3=set()
dict={}

#get all CDR3s
file=open(args.input)
reader=csv.reader(file,delimiter="\t")
next(reader,None)
for line in reader:
    cdr3=line[3]
    if cdr3[0] == "C" and cdr3[len(cdr3) - 1] == "F":
        set_cdr3.add(cdr3)
file.close()


for c in set_cdr3:
    dict[c]=0

#get number of reads

file=open(args.input)
reader=csv.reader(file,delimiter="\t")
next(reader,None)
for line in reader:
    nReads=int(line[2])
    cdr3=line[3]
    if cdr3[0]=="C" and cdr3[len(cdr3)-1]=="F":
        dict[cdr3]+=nReads
    else:
        print (cdr3)

file.close()

fileOut=open(args.out_prefix+".cdr3.csv","w")
fileOut.write("CDR3,nReads")
fileOut.write("\n")

for key, value in dict.items():
    fileOut.write(key+","+str(value))
    fileOut.write("\n")


fileOut.close()

# alpha diversity
fileOut=open(args.out_prefix+".alpha.diversity.csv","w")
fileOut.write(args.out_prefix+","+str(sdi(dict)))
fileOut.write("\n")



