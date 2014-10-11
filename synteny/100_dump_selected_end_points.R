#!/usr/bin/Rscript

#
# dump selected points to csv file
#

#spp='sorghum'

inpfile=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_selected_points2',spp)
outfile=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_selected_pts2.csv',spp)

load(inpfile)

if(file.exists(outfile)) {
    file.remove(outfile)
}

for (i in 1:(dim(x)[1]))
{
    rw = x[i,]
    write(paste(lapply(rw,as.character),collapse=','),file=outfile,append=T)
}
