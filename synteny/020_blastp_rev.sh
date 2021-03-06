#$ -S /bin/sh
#$ -N blastp
#$ -o ../scr/logs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ../scr/logs/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -cwd
#$ -l h_vmem=1G
###$ -q amd.q
#$ -t 1-20

#
# blastp model genome proteins against predicted oat proteins
#

taskid=`printf "%03d" $((SGE_TASK_ID-1))`
db='../../annotation_2013-11-05/scr/atlantica_annotation_20140226/tg7_proteins.fa'

#spp='brachy'
#query="../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/brachy_split/Bdistachyon_192_protein_primaryTranscriptOnly-${taskid}.fa"
#spp='rice'
#query="../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/rice_split/Osativa_204_protein_primaryTranscriptOnly-${taskid}.fa"
#spp='sorghum'
#query="../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/sorghum_split/Sbicolor_79_protein_primaryTranscriptOnly-${taskid}.fa"
spp='barley'
query="../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/barley_split/barley_mapped_prots-${taskid}.fa"

output="../scr/blastp_rbh_results/rev_${spp}-${taskid}.tsv"

blastp='/cm/shared/apps/BLAST/ncbi-blast-2.2.28+/bin/blastp'
outfmt='6'
ops='-evalue 1e-30'

${blastp} -db ${db} -query ${query} -out ${output} ${ops} -outfmt ${outfmt}
