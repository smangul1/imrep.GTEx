import Levenshtein
import argparse
import csv
import sys
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





def sumFreq(list,dict):
    s=0.0
    for l in list:
        if l in dict.keys():
            s+=dict[l]
    return s


ap = argparse.ArgumentParser()
ap.add_argument('input', help='CDR3 file from BCRSEQ')
ap.add_argument('output', help='mixcr  file with CDR3s')
args = ap.parse_args()




dict_counts_bcrseq={}
dict_counts_bcrseq_vdj={}


clonotypes_bcr=set()

vdj_bcr=set()




vdj_bcr=set()

with open(args.input) as csvFile:
    readCSV=csv.reader(csvFile,delimiter="\t")
    next(readCSV, None)
    for line in readCSV:

        cdr3=line[3]
        count=int(line[2])


        #if line[2]!="Out" and '*' not in line[1] and :
        if cdr3[0]=="C" and cdr3[len(cdr3)-1]=="W":
            clonotypes_bcr.add(cdr3)
            dict_counts_bcrseq[cdr3]=0
            if line[5]!="unresolved" and line[7]!="unresolved":
                V=line[5].replace('0','').split("-")[0]
                J=line[7].replace('0','').split("-")[0]
                vdj_bcr.add(V+"-"+J)
                dict_counts_bcrseq_vdj[V+"-"+J]=0


with open(args.input) as csvFile:
    readCSV=csv.reader(csvFile,delimiter="\t")
    next(readCSV, None)
    for line in readCSV:

        cdr3=line[3]
        count=int(line[2])

        if cdr3[0] == "C" and cdr3[len(cdr3) - 1] == "W":
            dict_counts_bcrseq[cdr3]+=count
            if line[5]!="unresolved" and line[7]!="unresolved":
                V=line[5].replace('0','').split("-")[0]
                J=line[7].replace('0','').split("-")[0]

                dict_counts_bcrseq_vdj[V+"-"+J]+=count









out_bcrseq=open(args.output,"w")

out_bcrseq.write("CDR3,nReads,FREQ\n")

for key,value in dict_counts_bcrseq.items():
    s=sum(dict_counts_bcrseq.values())
    freq=value/float(s)
    out_bcrseq.write(key+","+str(value)+","+str(freq))
    out_bcrseq.write("\n")

out_bcrseq.close()

