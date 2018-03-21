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
ap.add_argument('output', help='imrep  file with CDR3s')
args = ap.parse_args()



dict_counts_mixcr={}
dict_counts_mixcr_vdj={}



clonotypes_mixcr=set()


vdj_mixcr=set()



#mixcr -----------------------------------------------------------------------------
#0       2283    0.9921773142112125      TGTCAGCAGTATGGTAACTTCCCCCTCACTTTC       JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ       IGKV3-20*00(174.9)              IGKJ4*00(97.9)  IGKC*00(98.7)   324|343|368|0|19|SG340A|81.0            21|30|58|24|33||45.0
        #TGTCAGCAGTATGGTAACTTCCCCCTCACTTTC       41                                                               CQQYGNFPLTF             :::::::::0::-5:::::-1::33:::


with open(args.input) as csvFile:
    readCSV=csv.reader(csvFile,delimiter="\t")
    next(readCSV, None)
    for line in readCSV:
        if "IGH" in line[5]:
            cdr3=line[32]
            count=int(line[1])
            clonotypes_mixcr.add(cdr3)
            V=line[5].split("-")[0]
            J=line[7].split("*")[0]
            vdj_mixcr.add(V+"-"+J)

            dict_counts_mixcr[cdr3]=0
            dict_counts_mixcr_vdj[V+"-"+J]=0


with open(args.input) as csvFile:
    readCSV=csv.reader(csvFile,delimiter="\t")
    next(readCSV, None)
    for line in readCSV:
        if "IGH" in line[5]:
            cdr3=line[32]
            count = int(line[1])
            clonotypes_mixcr.add(cdr3)
            V=line[5].split("-")[0]
            J=line[7].split("*")[0]
            vdj_mixcr.add(V+"-"+J)

            dict_counts_mixcr[cdr3]+=count
            dict_counts_mixcr_vdj[V+"-"+J]+=count







out_mixcr=open(args.output,"w")
out_mixcr.write("CDR3,nReads,FREQ\n")

for key,value in dict_counts_mixcr.items():
    s=sum(dict_counts_mixcr.values())
    freq=value/float(s)
    out_mixcr.write(key+","+str(value)+","+str(freq))
    out_mixcr.write("\n")

out_mixcr.close()