while read line; do SRA=$(echo $line | awk '{print $1}'); new=$(echo $line | awk '{print $2}'); echo $SRA,$new; mv mixcr_${SRA}.txt mixcr_${new}.txt;done<../../sample_SRA_bio_new.txt

