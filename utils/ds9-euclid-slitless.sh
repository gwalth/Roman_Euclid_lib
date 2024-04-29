#!/bin/sh

# last modified 4/19/24

################################
# input: 
#
# e.g.
#    cd /Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11/Prep
#    ds9-euclid.sh 1
#
################################

opt="-lock frame wcs -scale mode zscale -match scale -single
-zoom to fit -scale limits -1 3"

ds9 \
 -mosaic Euclid_FRAME$1_DET??_slitless_final.fits \
 -frame 2 Euclid-total_ref.fits \
 $opt

