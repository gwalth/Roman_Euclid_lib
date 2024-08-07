import math, glob, os, sys

from astropy.table import Table
import astropy.io.fits as pyfits
import astropy.wcs as pywcs

import numpy as np
from scipy import integrate, interpolate
from scipy import linalg
from scipy.optimize import leastsq, least_squares

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from matplotlib.colors import Normalize


#import grizli_aXe_glw
from grizli_aXe_glw import (
     aXeConf_glw, 
     EuclidData, 
     EuclidFilters, 
     zodi_subtract, 
     insidePolygon,
     flam_to_fnu, 
     get_range,
)
from grizli_functions import euclid_wcs, remove_basic_wcs, stats

import matplotlib as mpl    
mpl.rcParams['font.size'] = 8



simdata_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/Prep/"
raw_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/Raw/EuclidSIMS/FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/frame1/"

####################################
####################################
####################################
#gain = 3500.0
#gain = 0.5
#gain = 1.0
gain = 6.1
#aper = 10
aper = 20.
#aper = 30.
#aper = 50.


det = '11'
frame = 1

# invesitgate the brightest
#Nfirst = 6
#Nfirst = 12
#Nfirst = 18
#Nfirst = 24
#Nfirst = 36
Nfirst = -1 # all


lw = 0.8
alpha = 0.6


verb = 1

gfit = False

cmap = cm.viridis
#cmap = cm.gist_ncar
#cmap = cm.jet

ymin = 1.
ymax = -1.

wmin = 12000
wmax = 18450

skip_indicies = [
]
####################################
####################################
####################################

#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv1/"
#pre_process_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS000p0/"
conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv4/"
pre_process_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv4/"

#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS000p0/"    # Frame1   
#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS180p4/"    # Frame2
#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS180p0/"    # Frame3
#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS000m4/"    # Frame4


#roman_path = prefix + "roman_setup/"
#simdata_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11/Prep/"
simdata_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/Prep/"

#conf = "config_axesim_det_11.conf"
#conf = "NISP_RGS000_12.conf"
#conf = "NISP_RGS000_34.conf"
#conf = "NISP_RGS180p4_12.conf"
#conf = "test_NISP_RGS180p4_12.conf"
#conf = "test2_NISP_RGS180p4_12.conf" # RGS180p4
#conf = "NISP_RGS180p4_12_final.conf" # RGS180p4
#conf = "test_NISP_RGS000m4_12.conf" # RGS000m4
#conf = "test3_NISP_RGS000_12.conf"
conf = "NISP_RGS000p0_12.conf"

sens_file = 'SENS_A.fits'
#input_spec_list = "input_spectra.lis"

direct_file = simdata_path + "Euclid-NISP_H_ref.fits"
#slitless_file = simdata_path + "NISP_SLITLESS_FRAME1_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits"   # RGS000p0
#slitless_file = simdata_path + "NISP_SLITLESS_FRAME2_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits"   # RGS180p4
#slitless_file = simdata_path + "NISP_SLITLESS_FRAME3_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits"   # RGS180p0
#slitless_file = simdata_path + "NISP_SLITLESS_FRAME4_FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11.fits"   # RGS000m4


slitless_file = simdata_path + "NISP_SLITLESS_FRAME%i_FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17.fits" % (frame)



#idx = 0
#xc = cat[idx]['X_IMAGE']
#yc = cat[idx]['Y_IMAGE']
#mag = cat[idx]['MAG_AUTO']
#spectemp = cat[idx]['SPECTEMP']

#print(xc, yc, mag)

##################################################

#aXe = aXeConf_glw(conf_file=conf_path + conf)

#print(aXe.conf)
#aXe.get_beam(xc, yc)
#aXe.load_files(direct_file, slitless_file, offset=0)
#aXe.load_files(direct_file, slitless_file, offset=-1024.)

#aXe.load_sensitivity(conf_path+sens_file)

#idx = 0
#xc = cat[idx]['X_IMAGE']
#yc = cat[idx]['Y_IMAGE']
#mag = cat[idx]['MAG_AUTO']
#spectemp = cat[idx]['SPECTEMP']

#print(xc, yc, mag)

##################################################



aXe = aXeConf_glw(conf_file=conf_path + conf)

print(aXe.conf)
#aXe.get_beam(xc, yc)
#aXe.load_files(direct_file, slitless_file, offset=0)
#aXe.load_files(direct_file, slitless_file, offset=-1024.)
aXe.load_files(direct_file, slitless_file, pre_process_path=pre_process_path, offset=1024, fix_mosa=1, det=det,
remove_dc=0, correct_gain=1, diag=1, gain=1.0)

aXe.load_sensitivity(conf_path+sens_file)


# In[47]:


euc_dat = EuclidData(
    catalog1 = "FSpatch_mod3_16183_TAcalib_newGRID_V1.fits",
    catalog2 = "NisInputConfiguration_%i_catalog_%s.cat" % (frame, det),
    model_spectra = "catalog_%s_%s.%i.spc.fits" % (det, det, frame),
    euclud_sims_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/Raw/EuclidSIMS/",
    catalog1_path = "FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/",
    catalog2_path = "FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/frame%i/Catalogs/" % (frame),
    model_spectra_path = "FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/frame%i/Input_Spectra/" % (frame),
)


# In[48]:


# slitless is 2048 x 2048
# reference is 10000 x 10000

#catalog_file = simdata_path + "MOT_SCA1_roll_0_dither_0x_0y_cut_zcut_cat14-18.fits"
#catalog_file = simdata_path + "Euclid_total_ref.cat"
catalog_file = simdata_path + "Euclid-NISP_H_ref.cat"

print(catalog_file)
#primer = Table.read(catalog_file, format='fits') 
cat = Table.read(catalog_file, format='ascii.sextractor') 

#print(cat.colnames)
#print(len(cat))


points = np.array([[a,d] for a,d in cat["X_WORLD","Y_WORLD"]])
print(points)

A = aXe.footprint_slit

#print(A)
#inside_det = []
#for p0 in points:
#    #print(p0)
#    inside = insidePolygon(A, p0)
#    #print(inside)
#
#    inside_det.append(inside)

inside_det = np.array([insidePolygon(A, p0) for p0 in points]) 
inside_det = np.multiply(inside_det, 1)

#########################################
# return indices of filter (boolean mask)
#########################################
det_indices = np.where(inside_det)[0]
#print(det_indices)
#cat_ids = [cat['NUMBER'][i] for i in det_indices]
#print(cat_ids)



#points_culled = points[inside_det]    
points_culled = points[det_indices]    
#print(points_culled)
#print(len(points_culled))
#skycoords = pw.wcs_pix2world(pixcoords,1)
pixcoords = aXe.wcs_slit.wcs_world2pix(points_culled,1)
#print(pixcoords)
#xc, yc = pixcoords[0]



#det_mag_idx = [[i,m] for i,m in enumerate(cat[det_indices]["MAG_AUTO"])]
#det_mag_idx_sorted = sorted(det_mag_idx, key=lambda obj: obj[1])
#det_idx = [idx for idx,mag in det_mag_idx_sorted]
#det_mag = [mag for idx,mag in det_mag_idx_sorted]

#print(det_mag_idx_sorted)
#print(det_idx)
#print(det_mag)
#print(len(det_mag_idx_sorted))




#print(points_culled[det_idx])
#print(pixcoords[det_idx])
#print(cat[det_indices]["MAG_AUTO"][det_idx])





#ratio_spectra = []



#subcat = cat["MAG_AUTO"][det_idx]
#print(subcat)
#print("I am here bitches!")


print(euc_dat.cat1)
print(euc_dat.cat1["Z_OBS"])
#zobs = euc_dat.cat1["Z_OBS"][det_indices][det_idx]
zobs = euc_dat.cat1["Z_OBS"][det_indices]
print(zobs)
print(euc_dat.cat1.keys())
#sys.exit()

#line = 6564.614
lines = [5008.239, 6564.614]


# https://classic.sdss.org/dr7/products/spectra/vacwavelength.php
all_lines = [
    # line      air       vacuum
    ["H-beta",  4861.363, 4862.721],
    ["[O III]", 4958.911, 4960.295],
    ["[O III]", 5006.843, 5008.239],
    ["[N II]",  6548.05,  6549.86],
    ["H-alpha", 6562.801, 6564.614],
    ["[N II]",  6583.45,  6585.27],
    ["[S II]",  6716.44,  6718.29],
    ["[S II]",  6730.82,  6732.68],
]


sci = aXe.slitless[0]
vmin, vmax = get_range(sci, sigma=5.) 

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.imshow(sci, origin="lower", vmin=vmin, vmax=vmax)



#for i,(p0, c0) in enumerate(zip(points_culled[det_idx][:Nfirst], pixcoords[det_idx][:Nfirst])):
for i,(p0, c0) in enumerate(zip(points_culled, pixcoords)):


    #euc_dat.get_spectra_ids()

    a0, d0 = p0
    xc, yc = c0

    print(xc, yc)

    sid, ind = euc_dat.find_source(a0, d0)
    print(sid,ind)
    model_wav, model_flux = euc_dat.get_model_spectra(a0, d0)
    #euc_dat.plot_model_spectra()



    #mag = det_mag[i]

    z = euc_dat.cat1["Z_OBS"][ind]
    #z = zobs[i]






    #aXe.get_beam_invert_cutout(xc, yc, aper, beam='A', offset=0.0)
    aXe.get_beam_invert_cutout(xc, yc, aper, beam='A', offset=0.0, dy0=-276.119403, dy1=320.895522)
    #aXe.extract_cutout(disp_axis=1, spatial_axis=0, flip_axis=1)
    #print(aXe.pix)
    #print(aXe.wav)
    #print(aXe.y_dist)

    w0 = np.min(aXe.wav)
    w1 = np.max(aXe.wav)

    interp_optical_model_disp = interpolate.interp1d(aXe.wav, aXe.pix)
    interp_optical_model_y_dist = interpolate.interp1d(aXe.wav, aXe.y_dist)


    for line in lines:

        wobs = line * (1 + z)
        print("wav(obs) = %.2f" % (wobs))

        if wobs > w0 and wobs < w1:

            y_pix_wav = interp_optical_model_disp(wobs)
            x_pix_wav = interp_optical_model_y_dist(wobs)
            print(x_pix_wav)
            print(y_pix_wav)

            ax1.scatter(x_pix_wav, y_pix_wav, marker="o", fc="None", ec="red", s=100)
            ax1.text(x_pix_wav, y_pix_wav, "%.1f" % (line), c="red")

    #flux = aXe.spectrum[2]
    #err = aXe.spectrum[3]
    #wav = aXe.wav
    #profile = aXe.profile


    #flux /= gain
    #err /= gain
    print()



plt.show()
