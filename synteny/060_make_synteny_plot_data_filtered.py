#!/usr/bin/python

'''
create plot data to
compare putative ortholog positions in oat and model genome
after filtering

run on bert
'''

import sys

#spp = 'brachy'
#spp = 'rice'
#spp = 'sorghum'
#spp = 'barley'
spp = sys.argv[1]

oatpos='../../zipper_003/scr/vcf_fixed/gene_snppav_posn_high.csv'
oatsize='../../zipper_003/scr/vcf_fixed/lg_sizes.csv'

spppos='../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/%s_gene_posn.csv'%spp
sppsize='../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/%s_chrm_sizes.csv'%spp

hits='../scr/synteny/%s_synteny_filtered055.tsv'%spp
out='../scr/synteny/%s_synteny_plot_filtered055.csv'%spp

#load chrm / lg sizes
f = open(oatsize)
oatlg = f.readline().strip().split(',')
f.close()
oatlg = [float(x) for x in oatlg]
oatcum = [sum(oatlg[:i]) for i in xrange(len(oatlg))]

f = open(sppsize)
spplg = f.readline().strip().split(',')
f.close()
spplg = [float(x) for x in spplg]
sppcum = [sum(spplg[:i]) for i in xrange(len(spplg))]

#load gene positions
oat = {}
f = open(oatpos)
for line in f:
    tok = line.strip().split(',')
    oat[tok[0]] = [int(tok[1]),float(tok[2])]
f.close()

sppp = {}
f = open(spppos)
for line in f:
    tok = line.strip().split(',')
    sppp[tok[0]] = [int(tok[1]),float(tok[2])]
f.close()

#process hits
f = open(hits)
fout = open(out,'wb')
for line in f:
    tok = line.strip().split('\t')
    oatid = tok[0]
    
    if spp == 'barley':
        sppid = tok[1] #barley
    else:
        sppid = tok[1].split('.')[0] #for rice, brachy, sorghum
    
    #not all oat genes have positions
    if not oatid in oat: continue
    
    if not sppid in sppp:
        #print sppid
        continue
    
    #calculate x,y position of the hit
    #spp pos => x
    #oat pos => y
    #oatgene = oat[oatid]
    #sppgene = sppp[sppid]
    #y = float(oatgene[1]) + float(oatcum[oatgene[0]-1])
    #x = float(sppgene[1]) + float(sppcum[sppgene[0]-1])
    
    oat_lg = oat[oatid][0]
    oat_cm = oat[oatid][1]
    spp_lg = sppp[sppid][0]
    spp_cm = sppp[sppid][1]
    
    y = float(oat_cm) + float(oatcum[oat_lg-1])
    x = float(spp_cm) + float(sppcum[spp_lg-1])
    
    row = [str(x),str(y),oatid,str(oat_lg),str(oat_cm),sppid,str(spp_lg),str(spp_cm)]
    
    fout.write(','.join(row) + '\n')
    
fout.close()    
f.close()






















