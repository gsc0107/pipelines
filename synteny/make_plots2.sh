#!/bin/sh

# run on bert

for spp in barley brachy rice sorghum
do
./200_plot_orthologs.R ${spp}\
    ../scr/synteny/${spp}_coords_nofilter.csv\
    ../scr/synteny/${spp}_nofilter.png
done
