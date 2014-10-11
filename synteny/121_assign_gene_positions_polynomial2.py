#!/usr/bin/python

'''
run this on bert

using the defined syntenic regions and the fitted polynomials
assign positions to oat genes based on positions of orthologs in sequenced species

output a single best guess position, provided all estimates are on same LG
an within 10cM
'''

import sys
sys.path.append('/ibers/ernie/home/rov/python_lib')
from rjv.fasta import *
from rjv.fileio import *

oat_genes = '../../annotation_2013-11-05/scr/atlantica_annotation_20140226/tg7_proteins_high.fa'
outfile = '../../zipper_003/scr/vcf_fixed/gene_synteny_posn.csv'

#spp in any order
spplist = ['sorghum','rice','brachy'] 
#spplist = ['sorghum','rice','barley'] 

verbose = False

#load all oat high conf gene ids
all_oat = {}
for fa in next_fasta(oat_genes):
    all_oat[fa['id']] = {}
    
for spp in spplist:
    polyfile    = '../scr/synteny/%s_poly_coeffs.csv'%spp
    synfile     = '../scr/synteny/%s_synteny_end_points'%spp
    spp_posfile = '../../annotation_2013-11-05/scr/phytozome_9.0_2013-12-03/%s_gene_posn.csv'%spp
    hits        = '../scr/blastp_rbh_results/%s_rbh.tsv'%spp
    selectfile  = '../scr/synteny/%s_selected_pts2.csv'%spp #syntenic hits
    sppoutfile  = '../../zipper_003/scr/vcf_fixed/gene_synteny_posn_%s.csv'%spp

    #index oat genes with orthologs
    hit_uids = {}
    f = open(hits)
    for line in f:
        tok = line.strip().split('\t')
        uid = tok[0]
        
        if 'PACid' in tok[1]:
            #phytozome gene ids
            sppid = tok[1].split('|')[0].split('.')[0]
        else:
            #barley gene ids
            sppid = tok[1]
        hit_uids[uid] = sppid
    f.close()
    
    #load gene positions in spp
    spp_pos = {}
    f = open(spp_posfile)
    for line in f:
        tok = line.strip().split(',')
        uid = tok[0]
        lg = int(tok[1])
        pos = float(tok[2])
        spp_pos[uid] = [lg,pos]
    f.close()
    
    #load the selected homologs
    selected = []
    f = open(selectfile)
    for line in f:
        tmprow = line.strip().split(',')
        row = tmprow[2:]
        #      oatid  oatlg       oatpos        sppid  spplg       spppos        plot_x           plot_y
        row = [row[0],int(row[1]),float(row[2]),row[3],int(row[4]),float(row[5]),float(tmprow[0]),float(tmprow[1])]
        selected.append(row)
    f.close()
    
    #load syntenic regions
    syn_regs = []
    f = open(synfile)
    while True:
        row1 = f.readline()
        if row1 == '': break
        
        row1 = row1.strip().split(',')[2:]
        #       oatid   oatlg        oatpos         sppid   spplg        spppos
        row1 = [row1[0],int(row1[1]),float(row1[2]),row1[3],int(row1[4]),float(row1[5])]
        row2 = f.readline().strip().split(',')[2:]
        row2 = [row2[0],int(row2[1]),float(row2[2]),row2[3],int(row2[4]),float(row2[5])]
        
        #make sure LGs are the same!
        assert row1[4] == row2[4]
        
        #make sure row1 is 5' end of region
        if row1[5] > row2[5]: row1,row2 = row2,row1
        
        syn_regs.append([row1,row2])
        
    f.close()
    
    #load polynomial coefficients
    poly_coeffs = []
    for tok in next_line(polyfile):
        poly_coeffs.append([float(x) for x in tok])
    
    fout = open(sppoutfile,'wb')
    
    #for each oat gene
    for uid in all_oat.iterkeys():
        #ignore if no spp homolog
        if not uid in hit_uids:
            if verbose: print 'no_hit'
            continue
        
        #find the spp homolog
        sppid = hit_uids[uid]

        #ignore if no position for the spp gene
        #(ie not on main chromosome scaffold)
        if not sppid in spp_pos:
            if verbose: print 'no_spp_pos'
            continue
        
        #find position in spp
        spplg = spp_pos[sppid][0]
        spppos = spp_pos[sppid][1]
        
        #find matching syntenic region(s)
        syn_offset = None
        for i,reg in enumerate(syn_regs):
            assert reg[0][4] == reg[1][4]
            if spplg != reg[0][4]: continue #wrong LG
            if spppos < reg[0][5] or spppos > reg[1][5]: continue #wrong region
            if syn_offset != None:
                #second match found, treat as no match
                #print syn_regs[i][0][4],syn_regs[i][0][5],syn_regs[i][1][5]
                #print syn_regs[syn_offset][0][4],syn_regs[syn_offset][0][5],syn_regs[syn_offset][1][5]
                #print spplg, spppos
                #print
                if verbose: print 'multiple_regions'
                syn_offset = None
                break
                
            #found one match, keep checking
            #a second match will invalidate this match
            syn_offset = i
            #if flag: print i
            
        #ignore if not in syntenic region
        #or if covered by two or more overlapping syntenic regions
        if syn_offset == None:
            if verbose: print 'zero_or_multiple_regions'
            continue
        
        #calculate the oat position using the polynomial
        new_lg = syn_regs[syn_offset][0][1]
        
        c = poly_coeffs[syn_offset]
        x = spppos
        new_pos = c[0] + c[1]*x + c[2]*x**2 + c[3]*x**3 

        all_oat[uid][spp] = [new_lg,new_pos]
        
        if verbose: print 'assigned_position'
        fout.write('%s,%d,%f\n'%(uid,new_lg,new_pos))
        
    fout.close()
        
fout = open(outfile,'wb')
for uid in all_oat.iterkeys():
    lg = [all_oat[uid][spp][0] for spp in all_oat[uid].keys()]
    pos = [all_oat[uid][spp][1] for spp in all_oat[uid].keys()]
    
    #print
    #print lg
    #print pos
    
    if len(lg) == 0: continue
    
    if len(lg) == 1:
        #only one spp available
        lg = lg[0]
        pos = pos[0]
    else:
        #must agree on LG
        if lg.count(lg[0]) != len(lg):
            continue
            
        lg = lg[0]
            
        #range must be < 10cM
        if abs(max(pos)-min(pos)) > 10.0:
            continue
            
        #pick median of 3 estimates
        if len(pos) == 3:
            pos.sort()
            pos = pos[1]
        else:
            #pick rice if available else sorghum
            if 'rice' in all_oat[uid]:
                pos = all_oat[uid]['rice'][1]
            else:
                pos = all_oat[uid]['sorghum'][1]
                
    #print lg,pos
                
    fout.write('%s,%d,%f\n'%(uid,lg,pos))
fout.close()
