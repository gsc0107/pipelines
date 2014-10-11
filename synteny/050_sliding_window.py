#!/usr/bin/python

'''
process gene positions using sliding window
to identify candidate syntenic regions
filter hits file to remove non-syntenic genes

this script filters based on density of putative orthologs
in the model genome
'''

import sys
import numpy as np

#spp = 'brachy'
#spp = 'rice'
#spp = 'sorghum'
#spp='barley'
spp = sys.argv[1]

oatpos='../../zipper_003/scr/vcf_fixed/gene_snppav_posn_high.csv'
oatsize='../../zipper_003/scr/vcf_fixed/lg_sizes.csv'

spppos='../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/%s_gene_posn.csv'%spp
sppsize='../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/%s_chrm_sizes.csv'%spp

hits='../scr/blastp_rbh_results/%s_rbh.tsv'%spp
out='../scr/synteny/%s_synteny_filtered050.tsv'%spp #filtered hits

outplot='../scr/synteny/%s_coords_nofilter.csv'%spp

bins = 50          #window size is genome_size/bins
samples = 500.0    #sample genes per window at this many placed per chromosome

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

#load RBH blast hits
f = open(hits)

noposn_oat = 0
noposn_spp = 0
bothposn = 0

fout = open(outplot,'wb')
data = []
for line in f:
    tok = line.strip().split('\t')
    oatid = tok[0]
    
    if spp == 'barley':
        sppid = tok[1] #barley
    else:
        sppid = tok[1].split('.')[0] #for rice, brachy, sorghum
    
    #not all oat genes have positions
    if not oatid in oat:
        #print oatid
        noposn_oat += 1
        continue
    
    if not sppid in sppp:
        #print sppid
        noposn_spp += 1
        continue
        
    bothposn += 1
    
    oat_lg = oat[oatid][0]
    oat_pos = oat[oatid][1]
    spp_lg = sppp[sppid][0]
    spp_pos = sppp[sppid][1]
    
    #      0     1     2      3       4      5
    row = [oatid,sppid,oat_lg,oat_pos,spp_lg,spp_pos]
    
    #row = [str(x) for x in row]
    data.append(row)
    
    y = float(oat_pos) + float(oatcum[oat_lg-1])
    x = float(spp_pos) + float(sppcum[spp_lg-1])
    
    fout.write(str(x) + ',' + str(y) + '\n')
    
f.close()
fout.close()

print 'oat genes with no map position',noposn_oat
print '%s genes with no genome position'%spp,noposn_spp
print 'genes with positions in both',bothposn

genome_size = sppcum[-1]
window_size = genome_size / float(bins) / 2.0 #half window size
print "window size:",window_size

#find the mean and stddev number of genes
#in the window size
cts = []
for i in xrange(len(oatlg)):
    for j in xrange(len(spplg)):
        #process one oat lg x model chrm at a time
        tmp_data = [x for x in data if x[2] == i+1 and x[4] == j+1]
        #print len(tmp_data)
        
        #scan the spp chromosome
        chrm_size = spplg[j]
        chrm_inc = chrm_size/samples
        chrm_posn = np.arange(0,chrm_size,chrm_inc)
        for x in chrm_posn:
            total = 0
            for y in tmp_data:
                if abs(x-y[5]) <= window_size:
                    total += 1

            cts.append(total)

mean_ct = sum(cts) / float(len(cts))
stddev_ct = np.std(cts)

print 'mean,stddev genes per bin',mean_ct,stddev_ct

fout = open(out,'wb')
#output only genes whose window contains an above average count
for i in xrange(len(oatlg)):
    for j in xrange(len(spplg)):
        tmp_data = [x for x in data if x[2] == i+1 and x[4] == j+1]
        #print len(tmp_data)
        
        for x in tmp_data:
            #count genes within +/- window_size of gene x
            #NB count includes self
            total = 0
            for y in tmp_data:
                if abs(x[5]-y[5]) <= window_size:
                    total += 1

            if float(total-1) > mean_ct: # subtract self from total
                fout.write('\t'.join([str(k) for k in x]) + '\n')
fout.close()
