#cd intermediate.files
#while read line; do SRA=$(echo $line | awk '{print $1}'); new=$(echo $line | awk '{print $2}'); echo $SRA,$new; mv mixcr_${SRA}.txt mixcr_${new}.txt;done<../../sample_SRA_bio_new.txt

ls *txt | awk -F ".txt" '{print $1}' >samples.txt

while read line
do
python ../../code/mixcr.extract.py ${line}.txt ${line}.clean.cdr3
done<samples.txt

