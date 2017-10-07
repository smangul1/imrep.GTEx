# Temporary notes

```
while read line; do ls /u/home/s/serghei/collab/GTEX/raw/afterQC_lostHuman_Fasta/* | grep ${line};done<spleen.txt  >spleen_samples.txt 
while read line; do ln -s $line data/;done<spleen_samples_10.txt
```

imrep
```
while read line; do echo "python ~/code/imrep/imrep.py --fastq ${line}.fasta.gz ${line}_o5.cdr3">run_${line}.sh;done<samples.txt
ls run*sh | awk '{i+=1;print "qsub -cwd -V -N imrep7"i" -l h_data=16G,time=06:00:00 "$1}' >all.sh

while read line; do echo "python ~/code/imrep/imrep.py -o 7 --fastq ${line}.fasta.gz ${line}_o5.cdr3">run_${line}.sh;done<samples.txt
while read line; do echo "python ~/code/imrep/imrep.py -o 9 --fastq ${line}.fasta.gz ${line}_o5.cdr3">run_${line}.sh;done<samples.txt
while read line; do echo "python ~/code/imrep/imrep.py -o 11 --fastq ${line}.fasta.gz ${line}_o5.cdr3">run_${line}.sh;done<samples.txt
while read line; do echo "python ~/code/imrep/imrep.py -o 11 --fastq ${line}.fasta.gz ${line}_o5.cdr3">run_${line}.sh;done<samples.txt
```

