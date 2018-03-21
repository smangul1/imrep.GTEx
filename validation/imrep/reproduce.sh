
#while read line; do SRA=$(echo $line | awk '{print $1}'); new=$(echo $line | awk '{print $2}'); echo $SRA,$new; mv ${SRA}.cdr3 imrep_${new}.cdr3;done<../../sample_SRA_bio_new.txt


ls *cdr3 | awk -F ".cdr3" '{print $1}'>samples.txt

while read line
do
python ~/collab/code/imrep/clonality.py ${line}.cdr3 ${line}
done<samples.txt
