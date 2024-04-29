#!/bin/sh

# last modified 4/19/24

################################
# input: 
#
# e.g.
#    cd /Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11/Prep
#    ds9-euclid.sh
#
################################

opt="-rgb lock colorbar yes -rgb lock scale yes -lock frame wcs -scale limits -0.1 2 -match scale -single"

# Euclid_FRAME1_DET11_slitless_final

ds9 \
 Euclid-VIS_ref.fits \
 Euclid-total_ref.fits \
 -rgb -blue  Euclid-NISP_Y_ref.fits \
      -green Euclid-NISP_J_ref.fits \
      -red   Euclid-NISP_H_ref.fits \
 $opt

