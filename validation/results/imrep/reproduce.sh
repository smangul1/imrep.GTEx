ls *cdr3 | awk -F "." '{print $1}' >samples.txt
while read line;do python ~/collab/code/imrep/clonality.py ${line}.cdr3 ${line};done<samples.txt 

