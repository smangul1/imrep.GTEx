mixcr=/u/home/n/ngcrawfo/project-zarlab/igor/imrep_revision/tools/mixcr/mixcr/mixcr

sample=$3
mkdir $sample
cd $sample


$mixcr align -p rna-seq -OallowPartialAlignments=true -r log_align.txt $1 $2 alignments.vdjca
$mixcr assemblePartial -f -p alignments.vdjca alignments_rescued_1.vdjca
$mixcr assemblePartial -f alignments_rescued_1.vdjca alignments_rescued_2.vdjca
$mixcr assemble -ObadQualityThreshold=0 -f -r log_assemble.txt alignments_rescued_2.vdjca mixcr.clns
$mixcr exportClones -f mixcr.clns mixcr.txt

