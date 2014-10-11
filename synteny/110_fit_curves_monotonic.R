#!/usr/bin/Rscript

# run on laptop with bert directories mounted
# fit monotonic 3rd degree polynomials to each syntenic region

library(MASS)
library('MonoPoly') #not available on bert
library('plotrix')

sc = 0.5
mr = 0.7
pngwidth=800
pngheight=800

args = commandArgs(trailingOnly = TRUE)

spp=args[1]

inpfile=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_selected_points2',spp)
load(inpfile)
sizes='/home/rov/rjv_files/bert_scratch/zipper_003/vcf_fixed/lg_sizes.csv'
chrm_sizes=sprintf('/home/rov/rjv_files/bert_scratch/annotation_2013-11-05/phytozome_9.0_2013-12-03/%s_chrm_sizes.csv',spp)

endpts=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_synteny_end_points',spp)
selectpts=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_selected_pts2.csv',spp)

outplot=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_fitted_curves.png',spp)

outfile=sprintf('/home/rov/rjv_files/bert_scratch/synteny_2014-07-08/synteny/%s_poly_coeffs.csv',spp)
if(file.exists(outfile)) {
    file.remove(outfile)
}

png(outplot,width=pngwidth,height=pngheight)

xlabel=sprintf('%s (bp)',spp)
ylabel='atlantica(cM)'

par(col='black')

endz = read.table(endpts,sep=',',header=F)
ptz = read.table(selectpts,sep=',',header=F)


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
#plot(x$V1,x$V2,xlab=xlabel,ylab=ylabel,xlim=c(0,vlines[length(vlines)]),ylim=c(0,hlines[length(hlines)]))
plot(ptz$V1,ptz$V2,xlab=xlabel,ylab=ylabel,xlim=c(0,vlines[length(vlines)]),ylim=c(0,hlines[length(hlines)]))
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


#ptz / endz
#1 xpos
#2 ypos
#3 oatid
#4 oatlg
#5 oatpos
#6 sppid
#7 spplg
#8 spppos

synx = vector()
syny = vector()


#fit a 3rd degree polynomial to each syntenic region
#save the coefficients to file in the same order as the regions
#are listed in the regions file
for( i in 1:dim(endz)[1])
{
    if( i %% 2 == 0) next

    i1 = i
    i2 = i1+1
    
    oatlg = endz$V4[i1]
    spplg = endz$V7[i1]
    sppmin = endz$V8[i1]
    sppmax = endz$V8[i2]
    #cat(i1,i2,oatlg,spplg,sppmin,sppmax,'\n')
    
    #get points matching oat and spp lgs
    subptz = ptz[ptz$V4 == oatlg & ptz$V7 == spplg,]
    #print(dim(subptz)[1])
    
    #get points within syntenic region in spp genome
    subptz = subptz[subptz$V8 >= sppmin & subptz$V8 <= sppmax,]
    #print(dim(subptz)[1])
    
    synx = append(synx,endz$V1[i1])
    syny = append(syny,endz$V2[i1])
    synx = append(synx,endz$V1[i2])
    syny = append(syny,endz$V2[i2])
    synx = append(synx,NA)
    syny = append(syny,NA)
    
    #if(dim(subptz)[1] <= 2) next
    
    #fit model directly to plotting coordinates and draw on the plot
    pfx = subptz$V1
    pfy = subptz$V2
    df = data.frame(pfx,pfy)
    df = df[with(df,order(pfx)),]
    pfx = df$pfx
    pfy = df$pfy
    pf = monpol(pfy~pfx,df)
    coeffs = c(as.numeric(coef(pf)[1]),
               as.numeric(coef(pf)[2]),
               as.numeric(coef(pf)[3]),
               as.numeric(coef(pf)[4]))
    yyy = evalPol(pfx, coeffs)
    lines(pfx, yyy,col=2,lwd=1)
    
    #fit model properly to chromosomal position (bp) versus map position (cM)
    pfx = subptz$V8
    pfy = subptz$V5
    df = data.frame(pfx,pfy)
    pf = monpol(pfy~pfx,df)

    coeffs = c(as.numeric(coef(pf)[1]),
               as.numeric(coef(pf)[2]),
               as.numeric(coef(pf)[3]),
               as.numeric(coef(pf)[4]))

    write(paste(lapply(coeffs,as.character),collapse=','),file=outfile,append=T)
    
    #lines(synx,syny,col=3)
    #print(oatlg)
    #print(spplg)
    #print(dim(subptz))
}

dev.off()
