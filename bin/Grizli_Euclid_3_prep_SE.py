#!/usr/bin/env python

import sys
import os
import yaml
import time
import pickle

import astropy
import astropy.io.fits as pyfits
#from astropy import wcs
from astropy.table import Table, unique, join, vstack
import astropy.units as u
from astropy.coordinates import SkyCoord



import numpy as np

import matplotlib as mpl    
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator

#from IPython.display import Image

mpl.rcParams['figure.figsize'] = (10.0, 6.0)
mpl.rcParams['font.size'] = 12
mpl.rcParams['savefig.dpi'] = 72

import grizli_functions
#from grizli_functions import wcs_pixel_scale, check_sims, create_circular_mask
from grizli_functions import (
        add_noise, 
        wcs_pixel_scale, 
        check_sims, 
        display_grizli,
        fake_euclid_direct,
        fake_euclid_ref,
        total_image,
        euclid_det, 
        read_slitless_headers,
        write_individual_slitless, 
        plot_slitless, 
        map_src_to_det,
)
print(grizli_functions.__file__)

T_START = time.time()




####################################
# paramters
####################################
# noise characteristics of "direct image"
# Scaramella et al. 2022 - Euclid preparation I. The Euclid Wide Survey
# RGS000, RGS180, RGS000_rot, RGS180_rot
spec_exptime = 574 # seconds
#spec_gain = None
spec_gain = None
spec_gain = 0.25
#spec_gain = 0.5
#spec_gain = 1.0
#spec_gain = 6.1
offset = 1024.
rot = 2 # 2nd rotation option
correct_dc = 0 # use dc file
fr_str = "FRAME"



id_key = "NUMBER"
mag_key = "MAG_AUTO"
flux_key = "FLUX_AUTO"
x_image_key = "X_IMAGE"
y_image_key = "Y_IMAGE"
x_world_key = "X_WORLD"
y_world_key = "Y_WORLD"
#id_key = "id"
#mag_key = "mag_auto"
#flux_key =  'nisp_h_tot_corr'
#x_image_key = "x_image"
#y_image_key = "y_image"
#x_world_key = "x_world"
#y_world_key = "y_world"


verb = 1
plot = 1
#yaml_file = "config.yaml"
yaml_file = sys.argv[1]
plot_labels = 1
plot_mag_cut = 16.5
zp_mag_cut = 22.0
search_rad = 0.6 # arcsec
####################################

#YAML_PATH = os.getcwd()
with open(yaml_file, 'r') as f:
    yaml_dict = yaml.safe_load(f)
    print(yaml_dict)

HOME_PATH = yaml_dict["HOME_PATH"]
print("HOME_PATH =", HOME_PATH)
root = yaml_dict["root"]
print("root =", root)
YAML_PATH = os.path.join(HOME_PATH, root)

#dir()
#all_slitless = yaml_dict["all_slitless"]
#all_zodi = yaml_dict["all_zodi"]
#spec_exptime = yaml_dict["spec_exptime"]
#spec_gain = yaml_dict["spec_gain"]
all_cat = yaml_dict["all_cat"]
all_ref_files = yaml_dict["all_ref_files"]
slitless_files = yaml_dict["slitless_files"]
catalog_files = yaml_dict["catalog_files"]
phot_mode = yaml_dict["phot_mode"] 


Euclid_bands = ['VIS','NIR_Y','NIR_J','NIR_H']
Euclid_bands_flux = ['TU_FNU_VIS_MAG', 'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 'TU_FNU_H_NISP_MAG'] 
# ## Add all of the header metadata needed for Grizli


#print(all_slitless)
#print(spec_exptime)



os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

suffix = "_final"

all_final_direct = []
all_final_slitless = []


# ## Subtract 1024 from science and sqrt(error), and adjust header keys

# In[ ]:

# find all of the frames in the slitless files, hopefully they are all the same
frames = []
for sf in slitless_files:
    L = sf.split("_")
    frames.append([L[i] for i,l in enumerate(L) if fr_str in l][0])

print(frames)




#for j,single_frame in enumerate(all_slitless):
#    for i,slitless_file in enumerate(single_frame):
#        print(slitless_file)
#sys.exit()

## Write individual files for each extension of the slitless spectra
# Grizli is easier to manage when writing out all of the files. 
# At some point we'll want to read the data extensions directly into Grizli, 
# this is currently a kludge.
#all_slitless = ["Euclid_FRAME%i" % (i+1) + "_DET%s_slitless.fits" for i,sf in enumerate(slitless_files)]
all_slitless = []
for i,sf in enumerate(slitless_files):

    file_str = "Euclid_%s" % (frames[i]) + "_DET%s_slitless.fits"

    tmp_all_slitless = write_individual_slitless(
        sf, 
        file_str = file_str,  
        rot=rot,
        spec_exptime = spec_exptime,
        spec_gain = spec_gain,
        correct_dc = correct_dc,
        offset = offset,
    )

    all_slitless.append(tmp_all_slitless)

print(all_slitless)
#all_zodi = [write_individual_slitless(sf, file_str ="Zodi_%s" % (frames[i]) + "_DET%s_slitless.fits", rot=rot) for i,sf in enumerate(zodi_files)]

os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

if plot:

    sf = all_slitless[0][0]
    pf = pyfits.open(sf)

    data = pf['SCI'].data
    #data = pf['SCI'].data - offset
    X = data.flatten()

    data = pf['ERR'].data
    Y = data.flatten()

    fig = plt.figure(figsize=(10,5))

    ax1 = fig.add_subplot(121)
    #ax1.hist(X, bins=20, range=(0,2000)) # Counts
    #ax1.hist(X, bins=200)
    ax1.hist(X, bins=100, range=(-1, 1))
    ax1.set_yscale("log")
    #ax1.set_xlabel("Counts")
    ax1.set_xlabel("ADU/s")
    ax1.set_title("SCI")
    #ax1.xaxis.set_major_locator(MultipleLocator(5))

    ax_ins = inset_axes(
        ax1, 
        width="100%", 
        height="100%",
        bbox_to_anchor=(0.125, 0.675, 0.3, 0.3), # right side
        bbox_transform=ax1.transAxes,
        borderpad=0
    )

    ax_ins.hist(X, bins=20)
    ax_ins.set_yscale("log")
    #ax_ins.set_xlabel("Counts")
    ax_ins.set_xlabel("ADU/s")

    ax_ins.xaxis.set_major_locator(MultipleLocator(10))
    #ax_ins.xaxis.set_minor_locator(MultipleLocator(0.1))
    #ax_ins.yaxis.set_major_locator(MultipleLocator(200))
    #ax_ins.yaxis.set_minor_locator(MultipleLocator(50))

    #ax_ins.tick_params(which='major', width=2)
    #ax_ins.tick_params(which='minor', width=1)


    ax2 = fig.add_subplot(122)
    #ax2.hist(Y,bins=20, range=(0,100)) # Counts
    ax2.hist(Y, bins=100)
    ax2.set_yscale("log")
    #ax2.set_xlabel("Counts")
    ax2.set_xlabel("ADU/s")
    ax2.set_title("ERR")

    plt.show()




#sys.exit()

###############
### SPECTRA ###
###############
for j,single_frame in enumerate(all_slitless):
    for i,slitless_file in enumerate(single_frame):
        print(slitless_file)
        
        new_slitless_file = slitless_file.replace(".fits",suffix+".fits")
        all_final_slitless.append(new_slitless_file)
            
        hdu = pyfits.open(slitless_file)
        hdu.info()
        
        #print(hdu[0].header)
        #print(hdu[1].header)
        
        print(hdu[0].header["GWA_POS"])
        print(hdu[0].header["GWA_ANG"])
        print(hdu[0].header["GWA_REF"])
        print(hdu[0].header["GWA_TILT"])
        print(hdu[1].header["DET_ID"])
        
        
        grism = hdu[0].header["GWA_POS"]
        det_id = hdu[1].header["DET_ID"]
        tilt = hdu[0].header["GWA_TILT"]
        
        #fake_filter = "-".join([grism,det_id])
        #print(fake_filter)

        if tilt < 0: ang = "m%i" % (tilt)
        elif tilt >= 0: ang = "p%i" % (tilt)

        fake_instr = "_".join(["NISP", grism, ang])
        print(fake_instr)
        
        #if all_zodi:
        #    zodi_file = all_zodi[i]
        #    print(zodi_file)
        #    hdu_zodi = pyfits.open(zodi_file)
        #    zodi = hdu_zodi[ext].data
        #else:
        #    zodi = 0.0
        
        ext = 0
        #hdu[ext].header['INSTRUME'] = 'NISP' 
        #hdu[ext].header['INSTRUME'] = 'NISP-GLWv1' 
        hdu[ext].header['INSTRUME'] = fake_instr
        # v2:
        # - optical model looks the same for each detector
        # - sensitivity for each detector (are they different?)
        hdu[ext].header['FILTER'] = 'RED'
        #hdu[ext].header['FILTER'] = fake_filter
        hdu[ext].header['EXPTIME'] = spec_exptime
        
        ext = 1
        #hdu[ext].header['INSTRUME'] = 'NISP'
        #hdu[ext].header['INSTRUME'] = 'NISP-GLWv1'
        hdu[ext].header['INSTRUME'] = fake_instr
        hdu[ext].header['FILTER'] = 'RED'
        #hdu[ext].header['FILTER'] = fake_filter
        hdu[ext].header['EXTVER'] = ext
        hdu[ext].header['EXPTIME'] = spec_exptime
        
        #sci = hdu[ext].data
        #hdu[ext].data = sci/spec_exptime/gain
        #hdu[ext].data = (sci-1024.)/spec_exptime/spec_gain
        #hdu[ext].data = (1.0*sci - 1.0*zodi)/spec_exptime/spec_gain
        
        #ext = 2
        #chi2 = hdu[ext].data
        #hdu[ext].data = np.sqrt(chi2)/spec_exptime/spec_gain
        
        hdu.writeto(new_slitless_file, overwrite=True, output_verify='fix')
        print("Writing",new_slitless_file)


# In[ ]:


print(all_final_direct)
print(all_final_slitless)


# ## Check the pixel scale of the images

# In[ ]:


os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

#for i in range(4):
#    for j in range(4):
#        num = '%i%i' % (i+1,j+1)
#        new_direct = "Euclid_DET%s%s.fits" % (num,suffix)
#        new_slitless = "Euclid_DET%s_slitless%s.fits"  % (num,suffix)

for direct in all_final_direct:        
    wcs_pixel_scale(direct)

for slitless in all_final_slitless:
    wcs_pixel_scale(slitless)


# In[ ]:


all_phot = []


for cat in all_cat:
    #print(cat)
    if phot_mode == "SExtractor": format="ascii.sextractor"
    elif phot_mode == "SEP": format = "fits"

    phot = Table.read(cat, format=format) # ref_cat in multimission
    all_phot.append(phot)
    print(cat," number_of_sources =",len(phot))


print()
print(phot.colnames)


# In[ ]:

if plot:
    fig = plt.figure(figsize=(10,6))
    
    k = 0
    #ax = []
    #if phot_mode == "SExtractor":

    if 1:


        for i,phot in enumerate(all_phot):
            
            cat = all_cat[i]
        
            #print(phot[mag_key])
            #print(phot[flux_key])
        
            if len(phot) > 0:    
                print(np.min(phot[mag_key]),np.max(phot[mag_key]))
                #print(np.min(phot[flux_key]),np.max(phot[flux_key]))
        
            
            ax = fig.add_subplot(2,3,k+1)
            ax.hist(phot[mag_key],range=(10,32),bins=44)
            #ax.hist(phot[mag_key],bins=20)
            #ax.hist(phot[flux_key],bins=20)
        
        
            ax.set_title(cat, fontsize=12)
            ax.set_xlabel("mag")
            ax.set_ylabel("N")
            ax.set_xlim(10,32)
            
            k += 1
    #elif phot_mode == "SEP":

        
    plt.show()
    
    
# In[ ]:


for phot in all_phot:
    print(phot)


# In[ ]:

if plot:

    
    fig = plt.figure(figsize=(10,6))
    
    k = 0
    for direct,phot in zip(all_ref_files,all_phot):
    
    
        hdu = pyfits.open(direct)
        img = hdu[1].data

        filt = phot[mag_key] < plot_mag_cut

        subphot = phot[filt]
        
        ax1 = fig.add_subplot(2,3,k+1)
        ax1.imshow(img,origin="lower", vmin=-0.1, vmax=0.1)
        ax1.scatter(subphot[x_image_key], subphot[y_image_key], fc="None", ec="r", s=50, lw=0.5, alpha=1.0)
    
        if plot_labels:
            for j,src in enumerate(subphot[mag_key]):
            #for j,src in enumerate(subphot[id_key]):
            #for j,src in enumerate(subphot[flux_key]):
                ax1.text(subphot[x_image_key][j], subphot[y_image_key][j], "%.2f" % (src), fontsize=8)

        ax1.set_title(direct, fontsize=12)
        
        k += 1
    plt.show()
    
    


# ## Read the source catalog
#os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

print(catalog_files)
primer = Table.read(catalog_files[0]) 
print(primer.colnames)
print(len(primer))


#filt = primer['RA'] < 200.
#print(filt)
#print(primer[filt])

#index = np.arange(len(primer))
#print(index[filt])
#primer.remove_row(index[filt][0])

# add magnitudes band into catalog
for bk,fk in zip(Euclid_bands,Euclid_bands_flux):
    fnu_Jy = primer[fk] # Jy    
    mab = -2.5*np.log10(fnu_Jy) + 8.90   # Jy --> AB mag            
    mab[np.isinf(mab)]=-99.
    primer[bk] = mab    



if plot:
    fig = plt.figure(figsize=(10,6))

all_phot_matched_clean = []
    
k = 0
for i,phot in enumerate(all_phot):
    
    #print(all_final_direct[i])
    print("Total number of sources found =",len(phot))
    print("Search primer and all_tbl for RA/DEC matches")
    
    if len(phot) > 0:
        #c_prime = SkyCoord(ra=primer["RA"], dec=primer["DEC"])
        #c_prime = SkyCoord(ra=primer["RA_MAG"], dec=primer["DEC_MAG"])
        c_prime = SkyCoord(ra=primer["RA_MAG"]*u.deg, dec=primer["DEC_MAG"]*u.deg)
        c_phot = SkyCoord(ra=phot[x_world_key], dec=phot[y_world_key])


        #idx, d2d, d3d = c_prime.match_to_catalog_sky(c_phot)
        idx, d2d, d3d = c_phot.match_to_catalog_sky(c_prime)

        filt = d2d < 1*u.arcsec

        if plot:
            ax1 = fig.add_subplot(2,3,k+1)
            ax1.hist(d2d[filt].value*3600.,bins=25)
            ax1.set_xlabel("Distance [arcsec]")
            ax1.set_title(all_ref_files[i], fontsize=12)
        
        #primer['idx'] = idx
        #primer['d2d'] = d2d

        phot['idx'] = idx
        phot['d2d'] = d2d
        #print(primer.colnames)
        
        #phot['idx'] = np.arange(len(phot))
        primer['idx'] = np.arange(len(primer))
        #print(all_tbl.colnames)
        #print(all_tbl['idx'])
        
        print("Join all_tbl and primer tables")
        match_tbl = join(phot, primer, keys='idx')
        #print(match_tbl)
        #print(match_tbl.colnames)
        
        print("Select only sources < %.2f arcsec" % (search_rad))
        #filt = match_tbl['d2d'] < 1.0*u.arcsec
        clean_filt = match_tbl['d2d'] < search_rad*u.arcsec
        match_clean_tbl = match_tbl[clean_filt]
        print(list(match_clean_tbl['idx']))
        print("Number of matches =",len(match_clean_tbl))
        print()
        
        if not i: print(match_clean_tbl.colnames)
        
    else:
        match_clean_tbl = []
    
    all_phot_matched_clean.append(match_clean_tbl)
    
    k+=1

if plot:
    plt.show()
        


# In[ ]:


if plot:
    
    fig = plt.figure(figsize=(10,6))
    
    k = 0    
    for direct,match_clean_tbl in zip(all_ref_files,all_phot_matched_clean):
    
    
        hdu = pyfits.open(direct)
        img = hdu[1].data
    
        ax1 = fig.add_subplot(2,3,k+1)
        ax1.imshow(img,origin="lower",vmin=-0.1,vmax=0.1)
        #ax1.scatter(phot[x_image_key],phot[y_image_key],fc="None",ec="r",s=60)

        if len(match_clean_tbl) > 0:

            filt = match_clean_tbl[mag_key] < plot_mag_cut
            
            sub_match_clean_tbl = match_clean_tbl[filt]

            ax1.scatter(
                sub_match_clean_tbl[x_image_key], 
                sub_match_clean_tbl[y_image_key], 
                lw=0.5,
                fc="None",
                ec="r",
                s=50
            )
        
            if plot_labels:
                for j,src in enumerate(sub_match_clean_tbl['SOURCE_ID']):
                    ax1.text(
                        sub_match_clean_tbl[x_image_key][j],
                        sub_match_clean_tbl[y_image_key][j], 
                        src,
                        fontsize=8
                    )
        ax1.set_title(direct, fontsize=12)

        k += 1
    plt.show()
    
    
print("N =",len(all_phot_matched_clean))
print()

print(all_phot_matched_clean)

verb = 0

for i,match_clean_tbl in enumerate(all_phot_matched_clean[:-1]):  
#for i,match_clean_tbl in enumerate(all_phot_matched_clean):  
    zps = []
    mags = []

    print(len(match_clean_tbl))
    
    if len(match_clean_tbl) > 0:
        print("i =",i)
        
        
        for k,row in enumerate(match_clean_tbl):
            
            #ind = np.argmax(match_clean_tbl[flux_key])

            #print("ind =",ind)
            dn = row[flux_key]
            mag = row[Euclid_bands[i]]
            flux = row[Euclid_bands_flux[i]]

            ## mag = -2.5 * log10(dn/exptime) + zp
            zp = mag + 2.5*np.log10(dn)
            
            
            mags.append(mag)
            zps.append(zp)
            if verb:
                print(k)
                print("dn =", dn)
                print("mag =", mag)
                print("flux =", flux)
                print("zp =", zp)
                print()

    zps = np.array(zps)
    mags = np.array(mags)

    print("Magnitudes:")
    print("min = %.3f" % (np.nanmin(mags)))
    print("max = %.3f" % (np.nanmax(mags)))
    print()

    print("Zeropoint:")
    print("mean   = %.3f" % (np.nanmean(zps)))
    print("median = %.3f" % (np.nanmedian(zps)))
    print("std    = %.3f" % (np.nanstd(zps)))
    print("min    = %.3f" % (np.nanmin(zps)))
    print("max    = %.3f" % (np.nanmax(zps)))

    filt = mags < zp_mag_cut

    print("Zeropoint (mag < %.2f):" % (zp_mag_cut)) 
    print("mean   = %.3f" % (np.nanmean(zps[filt])))
    print("median = %.3f" % (np.nanmedian(zps[filt])))
    print("std    = %.3f" % (np.nanstd(zps[filt])))
    print("min    = %.3f" % (np.nanmin(zps[filt])))
    print("max    = %.3f" % (np.nanmax(zps[filt])))


    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(mags, zps, s=0.5)
        ax.plot([zp_mag_cut, zp_mag_cut], [np.nanmin(zps), np.nanmax(zps)], "-", c="grey", label="Magnitude cut", zorder=10)
        ax.plot([np.nanmin(mags), np.nanmax(mags)], [np.nanmean(zps), np.nanmean(zps)], "--", c="tab:orange", label="Mean", zorder=10)
        ax.plot([np.nanmin(mags), zp_mag_cut], [np.nanmean(zps[filt]), np.nanmean(zps[filt])], "-", c="tab:orange", label="Mean (cut)", zorder=10)
        ax.plot([np.nanmin(mags), np.nanmax(mags)], [np.nanmedian(zps), np.nanmedian(zps)], "--", c="tab:green", label="Median", zorder=10)
        ax.plot([np.nanmin(mags), zp_mag_cut], [np.nanmedian(zps[filt]), np.nanmedian(zps[filt])], "-", c="tab:green", label="Median (cut)", zorder=10)
        ax.set_title(all_ref_files[i], fontsize=12)
        ax.legend()
        ax.set_xlabel("mag")
        ax.set_ylabel("zeropoint")

        plt.show()


# test what the rms predicts for the magnitude limit
rms = 0.15
zp = 25.13
mag = - 2.5*np.log10(5*rms) + zp
print(mag)


# ## Euclid object simulation
# [top](#Table-of-Contents)

# In[ ]:


#import importlib
#importlib.reload(grizli.model)


# In[ ]:


#import grizli.model


# In[ ]:


print(all_final_slitless)
print(all_ref_files)
#print(all_seg)
#print(all_phot_matched_clean)
print(all_phot_matched_clean[-1][id_key, x_image_key, y_image_key, mag_key])

T_END = time.time()

print()
print("Finished in %.1f seconds" % (T_END-T_START))
#sys.exit()

with open('all_phot_matched_clean.pickle', 'wb') as f:
    # Pickle the 'data' dictionary using the highest protocol available.
    pickle.dump(all_phot_matched_clean, f, pickle.HIGHEST_PROTOCOL)

yaml_dict['all_slitless'] = all_slitless
#yaml_dict['all_zodi'] = all_zodi
yaml_dict["all_final_slitless"] = all_final_slitless
yaml_dict['spec_exptime'] = spec_exptime
yaml_dict['spec_gain'] = spec_gain

os.chdir(YAML_PATH)
with open(os.path.basename(yaml_file), 'w',) as f:
    yaml.dump(yaml_dict, f, sort_keys=False)


