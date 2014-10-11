#!/usr/bin/python

'''
from the gff files generate files giving the gene positions
for brachy, rice, sorghum
'''

#inp = '../../annotation_2013-11-05/scr/scr/phytozome_9.0_2013-12-03/Bdistachyon_192_gene.gff3_fixchrm.gff3'
#out = '../../annotation_2013-11-05/scr/scr/phytozome_9.0_2013-12-03/brachy_gene_posn.csv'
#inp = '../../annotation_2013-11-05/scr/scr/phytozome_9.0_2013-12-03/Osativa_204_gene.gff3_fixchrm.gff3'
#out = '../../annotation_2013-11-05/scr/scr/phytozome_9.0_2013-12-03/rice_gene_posn.csv'
inp = '../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/Sbicolor_79_gene.gff3_fixchrm.gff3'
out = '../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/sorghum_gene_posn.csv'

f = open(inp)
fout = open(out,'wb')
for line in f:
    if line.startswith('#'): continue
    tok = line.strip().split('\t')
    if not len(tok) == 9: continue
    
    attr = tok[8].split(';')
    
    geneid = attr[0].split('=')[1]
    chrm = int(tok[0])
    pos = int(tok[3])
    
    outline = [geneid,str(chrm),str(pos)]
    
    fout.write(','.join(outline) + '\n')
    
fout.close()
f.close()
