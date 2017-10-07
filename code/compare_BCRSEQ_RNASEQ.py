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

def compare_cdr3(clonotypes1,clonotypes2):

    dList=[]
    
    cdr3List=[]

    for i in clonotypes1:
        t=set()
        t.clear()
        for j in clonotypes2:
            d=Levenshtein.ratio(i,j)

        
            if d>=0.9:
                t.add(d)
        if len(t)>0:
            dList.append(max(t))
            cdr3List.append(i)

    if len(dList)!=0:
        av=sum(dList)/len(dList)
    else:
        av=0

    return (len(dList), av,cdr3List)



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
dict_counts_imrep={}



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
        if line[2]!="Out" and '*' not in line[1]:
        #clonotypes1.append(line[0])
            clonotypes_bcr.add(line[1])
            
            dict_counts_bcrseq[line[1]]=0.0
            dict[line[1]]=0.0
            #print line[9],line[12],line[15], line
            if line[9]!="" and line[15]!="":
                V=line[9].replace('0','')
                J=line[15].replace('0','')
                vdj_bcr.add(V+"-"+J)
            



# freq is for nucl
with open(args.file1) as csvFile:
    readCSV=csv.reader(csvFile,delimiter="\t")
    next(readCSV, None)
    for line in readCSV:
        if line[2]!="Out" and '*' not in line[1]:
            #clonotypes1.append(line[0])
            dict[line[1]]+=float(line[7])
            dict_counts_bcrseq[line[1]]+=float(line[7])


#imrep
sumReadsImrep=0
with open(args.file2) as csvFile:
    readCSV=csv.reader(csvFile,delimiter="\t")
    next(readCSV, None)
    for line in readCSV:
    #if int(line[2])>2:
        clonotypes_imrep.add(line[0])
        dict_counts_imrep[line[0]]=int(line[2])
        sumReadsImrep+=int(line[2])
        if line[3].count(",")==0 and line[5].count(",")==0:
            
            
            V=line[3].split("/")[0]
            J=line[5]
            vdj_imrep.add(V+"-"+J)




with open(args.file3) as csvFile:
    readCSV=csv.reader(csvFile,delimiter="\t")
    next(readCSV, None)
    for line in readCSV:
        if "IGH" in line[5]:
            clonotypes_mixcr.add(line[32])
            V=line[5].split("-")[0]
            J=line[7].split("*")[0]
            vdj_mixcr.add(V+"-"+J)



bcr_vs_imrep=compare_cdr3(clonotypes_bcr,clonotypes_imrep)
imrep_vs_mixcr=compare_cdr3(clonotypes_imrep,clonotypes_mixcr)
bcr_vs_mixcr=compare_cdr3(clonotypes_bcr,clonotypes_mixcr)




freq_imrep=sumFreq (bcr_vs_imrep[2],dict)
freq_mixcr=sumFreq (bcr_vs_mixcr[2],dict)

max_freq_bcrseq=max(dict_counts_bcrseq.values())
max_freq_imrep=max(dict_counts_imrep.values())/float(sumReadsImrep)






vdj_bcr_imrep=len(vdj_imrep.intersection(vdj_bcr))
vdj_bcr_mixcr=len(vdj_mixcr.intersection(vdj_bcr))
vdj_imrep_mixcr=len(vdj_imrep.intersection(vdj_bcr))

sdi_bcrseq=sdi(dict_counts_bcrseq)
sdi_imrep=sdi(dict_counts_imrep)


print "args.file1, args.file2, args.file3,len(clonotypes_bcr), len(clonotypes_imrep),len(clonotypes_mixcr), bcr_vs_imrep[0],imrep_vs_mixcr[0],bcr_vs_mixcr[0],bcr_vs_imrep[1],imrep_vs_mixcr[1],bcr_vs_mixcr[1],len(vdj_bcr), len(vdj_imrep), len(vdj_mixcr),vdj_bcr_imrep,vdj_imrep_mixcr,vdj_bcr_mixcr, freq_imrep, freq_mixcr,sdi_bcrseq,sdi_imrep,max_freq_bcrseq,max_freq_imrep"

print args.file1, args.file2, args.file3,len(clonotypes_bcr), len(clonotypes_imrep),len(clonotypes_mixcr), bcr_vs_imrep[0],imrep_vs_mixcr[0],bcr_vs_mixcr[0],bcr_vs_imrep[1],imrep_vs_mixcr[1],bcr_vs_mixcr[1],len(vdj_bcr), len(vdj_imrep), len(vdj_mixcr),vdj_bcr_imrep,vdj_imrep_mixcr,vdj_bcr_mixcr, freq_imrep, freq_mixcr,sdi_bcrseq,sdi_imrep,max_freq_bcrseq,max_freq_imrep




