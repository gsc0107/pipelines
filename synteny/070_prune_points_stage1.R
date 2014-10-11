#!/usr/bin/Rscript

# this script is to manually prune all hits within a given square
# (ie chromosome-lg combination)
# left click on any dot with a square and all dots are removed
# right click when done
# right click immedaitely to remove none

# run on laptop with bert and bert_scratch dirs mounted
# launch R 
# set the variable spp
# eg spp='sorghum'
# then run this script using source('070_prune...')
# to maximise the plot window, run, maximise window, quit using right click
# then restart without closing the plot window

# if you make a mistake, simply right click then run the script again

library('plotrix')

sc = 0.5
mr = 0.5

#args = commandArgs(trailingOnly = TRUE)

#spp=args[1]
#spp='sorghum'

coords=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_synteny_plot_filtered055.csv',spp)
sizes='/home/rov/rjv_files/bert_scratch/zipper_003/vcf_fixed/lg_sizes.csv'
chrm_sizes=sprintf('/home/rov/rjv_files/bert_scratch/annotation_2013-11-05/phytozome_9.0_2013-12-03/%s_chrm_sizes.csv',spp)

outfile=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_selected_points',spp)

xlabel=spp
ylabel='atlantica(cM)'

par(col='black')
x = read.table(coords,sep=',')
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

boxed.labels(lg_xcoord1,lg_ycoord1,labels=lg_labels1)

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

boxed.labels(lg_xcoord2,lg_ycoord2,labels=lg_labels2)

while(T)
{
    #click on points with left button
    #click with right button to end
    ret = identify(x$V1,x$V2,plot=F,n=1)
    
    #print(ret)
    
    if(length(ret) == 0) break

    #get the oat and spp lg
    #delete all rows matching this combination
    oatlg = x[ret,'V4']
    spplg = x[ret,'V7']
    
    #print(oatlg)
    #print(spplg)

    x = x[!(x$V4==oatlg & x$V7==spplg),]
    
    #delete ret from the data frame
    #x = x[-ret,]
    
    #replot data
    #system(sprintf('echo %s >> %s',ret,outfile))
    plot(x$V1,x$V2,xlab=xlabel,ylab=ylabel,xlim=c(0,vlines[length(vlines)]),ylim=c(0,hlines[length(hlines)]))
    abline(h=hlines[1:length(hlines)-1])
    abline(v=vlines[1:length(vlines)-1])
    boxed.labels(lg_xcoord1,lg_ycoord1,labels=lg_labels1)
    boxed.labels(lg_xcoord2,lg_ycoord2,labels=lg_labels2)
}

save(x,file=outfile)
