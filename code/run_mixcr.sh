mixcr=/u/home/n/ngcrawfo/project-zarlab/igor/imrep_revision/tools/mixcr/mixcr/mixcr

sample=$3
mkdir $sample
cd $sample

#$mixcr align -p rna-seq -r log_align.txt $1 $2 alignments.vdjca
$mixcr assemblePartial -f -p alignments.vdjca alignments_rescued_1.vdjca
$mixcr assemblePartial -f alignments_rescued_1.vdjca alignments_rescued_2.vdjca
$mixcr assemble -f -r log_assemble.txt alignments_rescued_2.vdjca mixcr.clns
$mixcr exportClones -f mixcr.clns mixcr.txt
grep IGH mixcr.txt | grep -v "IGL" | grep -v "IGK" | cut -f 33 | grep -v "*" |  grep "^C" | grep "W$" | sort | uniq >mixcr_IGH_${sample}.cdr3
grep IGK mixcr.txt | grep -v "IGL" | grep -v "IGH" | cut -f 33 | grep -v "*" |  grep "^C" | grep "F$" | sort | uniq >mixcr_IGK_${sample}.cdr3 
grep IGL mixcr.txt | grep -v "IGK" | grep -v "IGH" | cut -f 33 | grep -v "*" |  grep "^C" | grep "F$" | sort | uniq >mixcr_IGL_${sample}.cdr3 

