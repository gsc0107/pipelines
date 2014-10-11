#!/cm/shared/apps/R/R-3.0.2/bin/Rscript

#
# make synteny plot of atlantica versus model species
#

library('plotrix')

sc = 0.5
mr = 0.7

args <- commandArgs(trailingOnly = TRUE)

spp       = args[1]
inpfile   = args[2]
outfile   = args[3]

filetype  = substr(outfile,nchar(outfile)-2,nchar(outfile))
pngwidth  = 800
pngheight = 800

coords = inpfile
#coords='../scr/synteny/tmp_synteny_plot.csv'
sizes='../../zipper_003/scr/vcf_fixed/lg_sizes.csv'
chrm_sizes=sprintf('../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/%s_chrm_sizes.csv',spp)

xlabel=sprintf('%s (bp)',spp)
ylabel='atlantica (cM)'

if(filetype == "pdf") {
    pdf(outfile)
} else {
    png(outfile,width=pngwidth,height=pngheight)
}

x = read.table(coords,sep=',')
hsizes=scan(sizes,sep=',')
vsizes=scan(chrm_sizes,sep=',')
hlines = cumsum(hsizes)
vlines = cumsum(vsizes)
par(mai=c(mr,mr,mr,mr))
par(cex=sc)
par(xaxs='i')
par(yaxs='i')
par(cex.lab=3.0)
par(cex.axis=3.0)
par(family="sans")
#par(cex.main=10.0)
#par(cex.sub=10.0)
plot(x$V1,x$V2,
     xlab=xlabel,ylab=ylabel,
     xlim=c(0,vlines[length(vlines)]),ylim=c(0,hlines[length(hlines)]))
abline(h=hlines[1:length(hlines)-1])
abline(v=vlines[1:length(vlines)-1])

#oat lg labels
n = length(hsizes)
lg_labels = rep(0,n)
lg_xcoord = rep(0,n)
lg_ycoord = rep(0.06,n)

for (i in 1:n)
{
    lg_labels[i] = sprintf('%d',i)
    lg_ycoord[i] = hlines[i] - hsizes[i]/2
    lg_xcoord[i] = 0
}

#boxed.labels(lg_xcoord,lg_ycoord,labels=lg_labels)

# lg labels
n = length(vsizes)
lg_labels = rep(0,n)
lg_xcoord = rep(0,n)
lg_ycoord = rep(0.06,n)

for (i in 1:n)
{
    lg_labels[i] = sprintf('%d',i)
    lg_xcoord[i] = vlines[i] - vsizes[i]/2
    lg_ycoord[i] = 0
}

#boxed.labels(lg_xcoord,lg_ycoord,labels=lg_labels)

dev.off()
