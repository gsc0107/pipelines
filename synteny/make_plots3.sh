#!/bin/sh

#run on laptop

for spp in barley brachy rice sorghum
do
./110_fit_curves_monotonic.R ${spp}
done
