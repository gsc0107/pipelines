#!/usr/bin/Rscript

# zoomed version of 080
# shows each lg/chrm combo in its own plot

# this script removes individual points from the synteny plot
# which might otherwise influence the curve fitting
# left click any points to remove them
# then right click to finish

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

oat_lgs = length(hsizes)
spp_chrms = length(vsizes)

x_filtered = data.frame()

par(xaxs='r')
par(yaxs='r')
par(mai=c(mr,mr,mr,mr))
par(cex=sc)

for(oatlg in 1:oat_lgs)
{
    for(sppchrm in 1:spp_chrms)
    {
        x_tmp = x[ x$V4 == oatlg & x$V7 == sppchrm, ]
        
        if(length(x_tmp) == 0) next
        
        vguides = seq(from=0.0,to=vsizes[sppchrm],length.out=100)
    
        plot(x_tmp$V8,x_tmp$V5,
             xlab=xlabel,ylab=ylabel,
             col='red',
             xlim=c(0,vsizes[sppchrm]),ylim=c(0,hsizes[oatlg]),
             panel.first=abline(v=vguides,col="gray"))

        while(T)
        {
            #click on points with left button
            #click with right button to end
            ret = identify(x_tmp$V8,x_tmp$V5,plot=F,n=1)
            
            if(length(ret) == 0) break

            #delete point
            x_tmp = x_tmp[-ret,]
            
            #replot data
            #system(sprintf('echo %s >> %s',ret,outfile))
            plot(x_tmp$V8,x_tmp$V5,
                 xlab=xlabel,ylab=ylabel,col='red',
                 xlim=c(0,vsizes[sppchrm]),ylim=c(0,hsizes[oatlg]),
                 panel.first=abline(v=vguides,col="gray"))
        }
        
        x_filtered = rbind(x_filtered,x_tmp)
    }
}

x = x_filtered
save(x,file=outfile)
