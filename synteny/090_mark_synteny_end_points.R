#!/usr/bin/Rscript

# run after 080... script
# from same R session

# mark individual syntenic regions
# mark the start and end of the region wrt model genome position
# can mark either end first

# script registered pairs of clicks, right click to end


library('plotrix')

sc = 0.5
mr = 0.5

#args = commandArgs(trailingOnly = TRUE)

#spp=args[1]
#spp='sorghum'

#coords=sprintf('%s_synteny_zipper.csv',spp)
sizes='/home/rov/rjv_files/bert_scratch/zipper_003/vcf_fixed/lg_sizes.csv'
chrm_sizes=sprintf('/home/rov/rjv_files/bert_scratch/annotation_2013-11-05/phytozome_9.0_2013-12-03/%s_chrm_sizes.csv',spp)

inpfile=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_selected_points2',spp)
outfile=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_synteny_end_points',spp)

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
if(file.exists(outfile)) {
    file.remove(outfile)
}

synx = vector()
syny = vector()
#fout = file(outfile)
while(T)
{
    #click on two points with left button
    #click with right button to end
    ret = identify(x$V1,x$V2,plot=F,n=2)
    
    if(length(ret) == 0) break

    #get the oat and spp plot position
    xpos1 = x[ret[1],'V1']
    ypos1 = x[ret[1],'V2']
    xpos2 = x[ret[2],'V1']
    ypos2 = x[ret[2],'V2']

    #note, start/end may be in either order 5',3' or 3',5'
    #depending on which order they were clicked on
    if(xpos1 <= xpos2)
    {
        row1 = x[ret[1],]
        row2 = x[ret[2],]
    }
    else
    {
        row1 = x[ret[2],]
        row2 = x[ret[1],]
    }
    
    #row1 = x[ret[1],]
    #row2 = x[ret[2],]
    write(paste(lapply(row1,as.character),collapse=','),file=outfile,append=T)
    write(paste(lapply(row2,as.character),collapse=','),file=outfile,append=T)
    #writeLines(paste(lapply(row1,as.character),collapse=','),fout)
    #writeLines(paste(lapply(row2,as.character),collapse=','),fout)
    #print(xpos1)
    #print(ypos1)
    #print(xpos2)
    #print(ypos2)
    
    synx = append(synx,xpos1)
    syny = append(syny,ypos1)
    synx = append(synx,xpos2)
    syny = append(syny,ypos2)
    synx = append(synx,NA)
    syny = append(syny,NA)
    
    #print(synx)
    #print(syny)
    
    #replot data
    #system(sprintf('echo %s >> %s',ret,outfile))
    par(col='black')
    plot(x$V1,x$V2,xlab=xlabel,ylab=ylabel,xlim=c(0,vlines[length(vlines)]),ylim=c(0,hlines[length(hlines)]))
    abline(h=hlines[1:length(hlines)-1])
    abline(v=vlines[1:length(vlines)-1])
    par(col='green')
    lines(synx,syny)
    #boxed.labels(lg_xcoord1,lg_ycoord1,labels=lg_labels1)
    #boxed.labels(lg_xcoord2,lg_ycoord2,labels=lg_labels2)
}

#close(fout)
