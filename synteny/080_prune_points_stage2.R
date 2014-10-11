#!/usr/bin/Rscript

# this script removes individual points from the synteny plot
# which might otherwise influence the curve fitting
# left click any points to remove them
# then right lcick to finish

# must be run after 070... script
# from the same R session

# if you make a mistake, right click then rerun the script

#this this from laptop with bert folders mounted, same as 070...R

library('plotrix')

sc = 0.5
mr = 0.5

#args = commandArgs(trailingOnly = TRUE)

#spp=args[1]
#spp='sorghum'

#coords=sprintf('%s_synteny_zipper.csv',spp)
sizes='/home/rov/rjv_files/bert_scratch/zipper_003/vcf_fixed/lg_sizes.csv'
chrm_sizes=sprintf('/home/rov/rjv_files/bert_scratch/annotation_2013-11-05/phytozome_9.0_2013-12-03/%s_chrm_sizes.csv',spp)

inpfile=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_selected_points',spp)
outfile=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_selected_points2',spp)


xlabel=spp
ylabel='atlantica(cM)'

par(col='black')
load(inpfile)
hsizes=scan(sizes,sep=',')
vsizes=scan(chrm_sizes,sep=',')
hlines = cumsum(hsizes)
vlines = cumsum(vsizes)
par(mai=c(mr,mr,mr,mr))
par(cex=sc)
par(xaxs='i')
par(yaxs='i')
plot(x$V1,x$V2,xlab=xlabel,ylab=ylabel,xlim=c(0,vlines[length(vlines)]),ylim=c(0,hlines[length(hlines)]))
abline(h=hlines[1:length(hlines)-1])
abline(v=vlines[1:length(vlines)-1])

#oat lg labels
n = length(hsizes)
lg_labels1 = rep(0,n)
lg_xcoord1 = rep(0,n)
lg_ycoord1 = rep(0.06,n)

for (i in 1:n)
{
    lg_labels1[i] = sprintf('%d',i)
    lg_ycoord1[i] = hlines[i] - hsizes[i]/2
    lg_xcoord1[i] = 0
}

#boxed.labels(lg_xcoord1,lg_ycoord1,labels=lg_labels1)

# lg labels
n = length(vsizes)
lg_labels2 = rep(0,n)
lg_xcoord2 = rep(0,n)
lg_ycoord2 = rep(0.06,n)

for (i in 1:n)
{
    lg_labels2[i] = sprintf('%d',i)
    lg_xcoord2[i] = vlines[i] - vsizes[i]/2
    lg_ycoord2[i] = 0
}

#boxed.labels(lg_xcoord2,lg_ycoord2,labels=lg_labels2)

while(T)
{
    #click on points with left button
    #click with right button to end
    ret = identify(x$V1,x$V2,plot=F,n=1)
    
    if(length(ret) == 0) break

    #get the oat and spp lg
    #delete point
    x = x[-ret,]
    
    #delete ret from the data frame
    #x = x[-ret,]
    
    #replot data
    #system(sprintf('echo %s >> %s',ret,outfile))
    plot(x$V1,x$V2,xlab=xlabel,ylab=ylabel,xlim=c(0,vlines[length(vlines)]),ylim=c(0,hlines[length(hlines)]))
    abline(h=hlines[1:length(hlines)-1])
    abline(v=vlines[1:length(vlines)-1])
    #boxed.labels(lg_xcoord1,lg_ycoord1,labels=lg_labels1)
    #boxed.labels(lg_xcoord2,lg_ycoord2,labels=lg_labels2)
}

save(x,file=outfile)
