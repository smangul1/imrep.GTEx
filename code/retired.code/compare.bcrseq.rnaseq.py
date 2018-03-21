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
ap.add_argument('file1', help='CDR3 file from BCRSEQ')
ap.add_argument('file2', help='imrep  file with CDR3s')
ap.add_argument('file3', help='mixcr  file with CDR3s')
args = ap.parse_args()




dict_counts_bcrseq={}
dict_counts_bcrseq_vdj={}

dict_counts_imrep={}
dict_counts_imrep_vdj={}


dict_counts_mixcr={}
dict_counts_mixcr_vdj={}


clonotypes_bcr=set()
clonotypes_imrep=set()
clonotypes_mixcr=set()

vdj_bcr=set()
vdj_imrep=set()
vdj_mixcr=set()

dict={}


vdj_bcr=set()

with open(args.file1) as csvFile:
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


with open(args.file1) as csvFile:
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









#imrep -----------------------------------------------------------------------------
#CGAGTPPASSCPSLGRGW      IGH     2       IGHV3   NA      IGHJ2

sumReadsImrep=0
with open(args.file2) as csvFile:
    readCSV=csv.reader(csvFile)
    next(readCSV, None)
    for line in readCSV:
        if line[1]=="IGH":
            cdr3=line[0]
            count=int(line[2])

            clonotypes_imrep.add(cdr3)
            dict_counts_imrep[cdr3]=count
            if line[3].count(";")==0 and line[5].count(";")==0: # not ambiguous
                V=line[3].split("/")[0]
                J=line[5]
                vdj_imrep.add(V+"-"+J)
                dict_counts_imrep_vdj[V+"-"+J]=0

with open(args.file2) as csvFile:
    readCSV = csv.reader(csvFile)
    next(readCSV, None)
    for line in readCSV:
        if line[1] == "IGH":
            cdr3 = line[0]
            count = int(line[2])
            if line[3].count(";") == 0 and line[5].count(";") == 0:  # not ambiguous
                V = line[3].split("/")[0]
                J = line[5]
                dict_counts_imrep_vdj[V + "-" + J] +=count



#mixcr -----------------------------------------------------------------------------
#0       2283    0.9921773142112125      TGTCAGCAGTATGGTAACTTCCCCCTCACTTTC       JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ       IGKV3-20*00(174.9)              IGKJ4*00(97.9)  IGKC*00(98.7)   324|343|368|0|19|SG340A|81.0            21|30|58|24|33||45.0
        #TGTCAGCAGTATGGTAACTTCCCCCTCACTTTC       41                                                               CQQYGNFPLTF             :::::::::0::-5:::::-1::33:::


with open(args.file3) as csvFile:
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


with open(args.file3) as csvFile:
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








########################################################################################

print ("bcr-imrep",len(clonotypes_bcr),len(clonotypes_imrep), len(clonotypes_bcr.intersection(clonotypes_imrep)))
print ("bcr-mixcr",len(clonotypes_bcr),len(clonotypes_mixcr), len(clonotypes_bcr.intersection(clonotypes_mixcr)))
print ("mixcr-imrep",len(clonotypes_mixcr),len(clonotypes_imrep), len(clonotypes_mixcr.intersection(clonotypes_imrep)))



out_bcrseq=open("bcrseqq.cdr3","w")

out_bcrseq.write("CDR3,nReads,FREQ\n")

for key,value in dict_counts_bcrseq.items():
    s=sum(dict_counts_bcrseq.values())
    freq=value/float(s)
    out_bcrseq.write(key+","+str(value)+","+str(freq))
    out_bcrseq.write("\n")

out_bcrseq.close()


out_imrep=open("imrep.cdr3","w")
out_imrep.write("CDR3,nReads,FREQ\n")

for key,value in dict_counts_imrep.items():
    s=sum(dict_counts_imrep.values())
    freq=value/float(s)
    out_imrep.write(key+","+str(value)+","+str(freq))
    out_imrep.write("\n")

out_imrep.close()

out_mixcr=open("mixcr.cdr3","w")
out_mixcr.write("CDR3,nReads,FREQ\n")

for key,value in dict_counts_mixcr.items():
    s=sum(dict_counts_mixcr.values())
    freq=value/float(s)
    out_mixcr.write(key+","+str(value)+","+str(freq))
    out_mixcr.write("\n")

out_mixcr.close()