#!/usr/bin/Rscript

# zoomed version of 090

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

oat_lgs = length(hsizes)
spp_chrms = length(vsizes)

if(file.exists(outfile)) {
    file.remove(outfile)
}

par(mai=c(mr,mr,mr,mr))
par(cex=sc)
par(xaxs='r')
par(yaxs='r')

for(oatlg in 1:oat_lgs)
{
    for(sppchrm in 1:spp_chrms)
    {
        x_tmp = x[ x$V4 == oatlg & x$V7 == sppchrm, ]
        
        if(length(x_tmp) == 0) next
        
        vguides = seq(from=0.0,to=vsizes[sppchrm],length.out=100)

        plot(x_tmp$V8,x_tmp$V5,
             col='red',
             xlab=xlabel,ylab=ylabel,
             xlim=c(0,vsizes[sppchrm]),ylim=c(0,hsizes[oatlg]),
             panel.first=abline(v=vguides,col="gray"))

        synx = vector()
        syny = vector()

        while(T)
        {
            #click on two points with left button
            #click with right button to end
            ret = identify(x_tmp$V8,x_tmp$V5,plot=F,n=2)
            
            if(length(ret) == 0) break

            #get the oat and spp plot position
            xpos1 = x_tmp[ret[1],'V8']
            ypos1 = x_tmp[ret[1],'V5']
            xpos2 = x_tmp[ret[2],'V8']
            ypos2 = x_tmp[ret[2],'V5']

            #note, start/end may be in either order 5',3' or 3',5'
            #depending on which order they were clicked on
            if(xpos1 <= xpos2)
            {
                row1 = x_tmp[ret[1],]
                row2 = x_tmp[ret[2],]
            }
            else
            {
                row1 = x_tmp[ret[2],]
                row2 = x_tmp[ret[1],]
            }
            
            write(paste(lapply(row1,as.character),collapse=','),file=outfile,append=T)
            write(paste(lapply(row2,as.character),collapse=','),file=outfile,append=T)
            
            synx = append(synx,xpos1)
            syny = append(syny,ypos1)
            synx = append(synx,xpos2)
            syny = append(syny,ypos2)
            synx = append(synx,NA)
            syny = append(syny,NA)
            
            #replot data
            par(col='black')
            plot(x_tmp$V8,x_tmp$V5,
                 xlab=xlabel,ylab=ylabel,col='red',
                 xlim=c(0,vsizes[sppchrm]),ylim=c(0,hsizes[oatlg]),
                 panel.first=abline(v=vguides,col="gray"))
                 
            par(col='green')
            lines(synx,syny)
            abline(v=synx)
            par(col='black')
        }
    }
}
