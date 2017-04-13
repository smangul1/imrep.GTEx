
import sys
from os import listdir, system
from Bio import SeqIO
from random import random, randint, choice
from numpy import cumsum
import numpy.random as npran


geneName = sys.argv[1]
nrGenerate = sys.argv[2]
tmpDir = sys.argv[3]
resultFile = sys.argv[4]
geomp = sys.argv[5]



rootDir = "IMGT/"
geneName2 = geneName
if geneName == "TRA":
    geneName2 = "TCRA"
if geneName == "TRB":
    geneName2 = "TCRB"

files = listdir(rootDir + geneName)

V = []
D = []
J = []
C = []


V_names = []
D_names = []
J_names = []
C_names = []


def populate_map(genemap, genenames, filename):
    handle = open(filename, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        aidi = record.id
        aidi = aidi.split("|")[1]
        aidi = aidi.replace("/", "")
        aidi = aidi.split("*")[0]
        genenames.append(aidi)
        genemap.append(record.seq.tostring().upper())
    handle.close()



for f in files:
    if f.startswith("vi"):
        populate_map(V, V_names, rootDir + geneName + "/" + f)
    elif f.startswith("di"):
        populate_map(D, D_names, rootDir + geneName + "/" + f)
    elif f.startswith("jay"):
        populate_map(J, J_names, rootDir + geneName + "/" + f)
    elif f.startswith("ci"):
        populate_map(C, C_names, rootDir + geneName + "/" + f)



genes = [V]
genenames = [V_names]
if D:
    genes.append(D)
    genenames.append(D_names)
if J:
    genes.append(J)
    genenames.append(J_names)
if C:
    genes.append(C)
    genenames.append(C_names)






sample = [0]
while 0 in sample:
    sample = npran.geometric(p=float(geomp), size=int(nrGenerate))



alltran = []

vdj_comb = set()
loginfo = []

for s in sample:
    i = 0
    regs = []
    regnames = []
    ready = False
    while not ready:
        vdj = []
        for j in range(len(genes)):
            pos = randint(0, len(genes[j]) - 1)
            vdj.append(pos)
        if tuple(vdj) in vdj_comb:
            continue
        if V:
            regs.append(genes[i][vdj[i]])
            regnames.append(genenames[i][vdj[i]])
            i += 1
        if D:
            regs.append(genes[i][vdj[i]])
            regnames.append(genenames[i][vdj[i]])
            i += 1
        if J:
            regs.append(genes[i][vdj[i]])
            regnames.append(genenames[i][vdj[i]])
            i += 1
        if C:
            regs.append(genes[i][vdj[i]])
            regnames.append(genenames[i][vdj[i]])
        valid_cdr3 = False
        z = 0
        while z < 10 and not valid_cdr3:
            nbregs = []
            for j in range(len(vdj) - 1):
                r = randint(0, 15)
                nbregs.append("".join([choice(["A", "C", "T", "G"]) for p in range(r)]))
            trans = regs[0]
            for x, y in zip(nbregs, regs[1:]):
                trans += x
                trans += y
            with open("%s/ref.fq" % tmpDir, "w") as f:
                f.write("@1\n")
                f.write("%s\n" % trans)
                f.write("+\n")
                f.write("%s\n" % (len(trans) * "^"))
            command = "java -jar LymAnalyzer_cmd_1.2.2.jar %s/ref.fq %s %s hs igor No No" % (tmpDir, tmpDir, geneName2)
            system(command)
            with open("%s/igor/igor_cdr3statistics.txt" % tmpDir) as f:
                a = f.readlines()
            a = map(lambda x: x.strip().split(), a)
            if len(a) == 2 and a[1][1] != "out_of_frame":
                valid_cdr3 = True
            z += 1
        if valid_cdr3:
            vdj_comb.add(tuple(vdj))
            vdjline = regnames
            vdjline.append(s)
            vdjline.append(s * 1.0 / sum(sample))
            vdjline.append(a[1][1])
            vdjline.append(a[1][2])
            loginfo.append(vdjline)
            ready = True
        else:
            i = 0
            regs = []
            regnames = []
    print "APPENDING TRANSQ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", trans
    for m in range(s):
        alltran.append(trans)

with open(resultFile, "w") as f:
    counter = 1
    for i in range(len(alltran)):
        f.write(">" + str(counter) + "\n")
        f.write(alltran[i] + "\n")
        counter += 1


with open(resultFile + ".info", "w") as f:
    for linfo in loginfo:
        f.write("%s\n" % "\t".join(map(str, linfo)))

