#$ -S /bin/sh
#$ -N blastp
#$ -o ../scr/logs/$JOB_NAME.$JOB_ID.$TASK_ID.out
#$ -e ../scr/logs/$JOB_NAME.$JOB_ID.$TASK_ID.err
#$ -cwd
#$ -l h_vmem=1G
#$ -l h_rt=3600
###$ -q amd.q
#$ -t 1-20

#
# blastp predicted oat proteins against rice or sorghum or brachy proteins
#

taskid=`printf "%03d" $((SGE_TASK_ID-1))`
query="../../annotation_2013-11-05/scr/atlantica_annotation_20140226/prot_split20/tg7_proteins-${taskid}.fa"

spp='brachy'
db='../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/Bdistachyon_192_protein_primaryTranscriptOnly.fa'
#spp='rice'
#db='../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/Osativa_204_protein_primaryTranscriptOnly.fa'
#spp='sorghum'
#db='../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/Sbicolor_79_protein_primaryTranscriptOnly.fa'
#spp='barley'
#db='../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/barley_mapped_prots.fa'

output="../scr/blastp_rbh_results/tmp_fwd_${spp}-${taskid}.tsv"

blastp='/cm/shared/apps/BLAST/ncbi-blast-2.2.28+/bin/blastp'
outfmt='6'
ops='-evalue 1e-30'

${blastp} -db ${db} -query ${query} -out ${output} ${ops} -outfmt ${outfmt}
