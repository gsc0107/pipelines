#!/usr/bin/python

'''
find reciprocal best hits
using forward and reverse blast searches

before running this, merge blast results into single files:

cat fwd_${spp}$-???.tsv > ${spp}_fwd.tsv
cat rev_${spp}$-???.tsv > ${spp}_rev.tsv
'''

import sys

#spp = 'barley'
#spp = 'brachy'
#spp = 'rice'
#spp = 'sorghum'

spp = sys.argv[1]

inpfwd = '../scr/blastp_rbh_results/%s_fwd.tsv'%spp
inprev = '../scr/blastp_rbh_results/%s_rev.tsv'%spp
out    = '../scr/blastp_rbh_results/%s_rbh.tsv'%spp
maxhits = 1

#process all fwd blastp hits, retain best hits for each oatid
fwd_map = {}
f = open(inpfwd)
for line in f:
    tok = line.strip().split('\t')
    oatid = tok[0]
    sppid = tok[1]
    evalue = float(tok[10])
    
    #update forward best hits
    if not oatid in fwd_map:
        fwd_map[oatid] = [evalue,[sppid]]
    elif evalue < fwd_map[oatid][0]:
        fwd_map[oatid] = [evalue,[sppid]]
    elif evalue == fwd_map[oatid][0]:
        fwd_map[oatid][1].append(sppid)
f.close()

#process all rev blastp hits, retain best hits for each sppid
rev_map = {}
f = open(inprev)
for line in f:
    tok = line.strip().split('\t')
    sppid = tok[0]
    oatid = tok[1]
    evalue = float(tok[10])
    
    #update reverse best hits
    if not sppid in rev_map:
        rev_map[sppid] = [evalue,[oatid]]
    elif evalue < rev_map[sppid][0]:
        rev_map[sppid] = [evalue,[oatid]]
    elif evalue == rev_map[sppid][0]:
        rev_map[sppid][1].append(oatid)
f.close()

#remove proteins with more than 'maxhits' hits
for oatid in fwd_map.iterkeys():
    if len(fwd_map[oatid][1]) > maxhits:
        fwd_map[oatid][1] = []

for sppid in rev_map.iterkeys():
    if len(rev_map[sppid][1]) > maxhits:
        rev_map[sppid][1] = []
        
#prune hits, retain only reciprocal best
fout = open(out,'wb')
for oatid in fwd_map.iterkeys():
    for sppid in fwd_map[oatid][1]:
        if sppid in rev_map:
            if oatid in rev_map[sppid][1]:
                fout.write('%s\t%s\t%e\n'%(oatid,sppid,fwd_map[oatid][0]))
fout.close()
