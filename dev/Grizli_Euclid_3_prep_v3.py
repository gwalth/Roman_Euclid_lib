import sys
import os
import yaml
import time
import pickle
import glob

import astropy
import astropy.io.fits as pyfits
#from astropy import wcs
from astropy.table import Table, unique, join, vstack
import astropy.units as u
from astropy.coordinates import SkyCoord

import numpy as np

import matplotlib as mpl    
import matplotlib.pyplot as plt
#from matplotlib.gridspec import GridSpec
#from matplotlib.ticker import MultipleLocator

#from IPython.display import Image

import grizli.prep

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

def add_grizli_headers(hdu, hdict, exts=[0,1,2,3]):
    
    for ext in exts:
        for key in hdict:
            hdu[ext].header[key] = hdict[key]
            
    return hdu




####################################
# paramters
####################################
#id_key = "NUMBER"
#mag_key = "MAG_AUTO"
#flux_key = "FLUX_AUTO"
#x_image_key = "X_IMAGE"
#y_image_key = "Y_IMAGE"
#x_world_key = "X_WORLD"
#y_world_key = "Y_WORLD"
id_key = "id"
mag_key = "mag_auto"
flux_key =  'nisp_h_tot_corr'
x_image_key = "x_image"
y_image_key = "y_image"
x_world_key = "x_world"
y_world_key = "y_world"


grism = "red"
prefix = "Euclid"
aper = 6
verb = 1
plot = 1
#yaml_file = "config.yaml"
yaml_file = sys.argv[1]
plot_labels = 1
#plot_mag_cut = 16.5
plot_mag_cut = 18.0
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
all_slitless = yaml_dict["all_slitless"]
all_zodi = yaml_dict["all_zodi"]
spec_exptime = yaml_dict["spec_exptime"]
spec_gain = yaml_dict["spec_gain"]
all_cat = yaml_dict["all_cat"]
all_ref_files = yaml_dict["all_ref_files"]
catalog_files = yaml_dict["catalog_files"]
phot_mode = yaml_dict["phot_mode"] 


Euclid_bands = ['VIS','NIR_Y','NIR_J','NIR_H']
Euclid_bands_flux = ['TU_FNU_VIS_MAG', 'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 'TU_FNU_H_NISP_MAG'] 
# ## Add all of the header metadata needed for Grizli


print(all_slitless)
print(spec_exptime)



os.chdir(os.path.join(HOME_PATH, root, 'Prep'))

suffix = "_final"

all_final_direct = []
all_final_slitless = []


# ## Subtract 1024 from science and sqrt(error), and adjust header keys

# In[ ]:


###############
### SPECTRA ###
###############
for slitless in all_slitless[0]:
    print(slitless)

    new_slitless = slitless.replace(".fits",suffix+".fits")
    all_final_slitless.append(new_slitless)



if 0:
    for slitless in all_slitless[0]:
        print(slitless)
    
        new_slitless = slitless.replace(".fits",suffix+".fits")
        all_final_slitless.append(new_slitless)
            
        hdu = pyfits.open(slitless)
        hdu.info()
    
        ext = 0
        hdu[ext].header['INSTRUME'] = 'NISP' 
        #hdu[ext].header['INSTRUME'] = 'NISP-GLWv3' 
        # v2:
        # - optical model looks the same for each detector
        # - sensitivity for each detector (are they different?)
        hdu[ext].header['FILTER'] = 'RGS000p0'
        hdu[ext].header['EXPTIME'] = spec_exptime
    
        ext = 1
        hdu[ext].header['INSTRUME'] = 'NISP'
        #hdu[ext].header['INSTRUME'] = 'NISP-GLWv3'
        hdu[ext].header['FILTER'] = 'RGS000p0'
        hdu[ext].header['EXTVER'] = ext
        hdu[ext].header['EXPTIME'] = spec_exptime
    
        sci = hdu[ext].data
        #hdu[ext].data = sci/spec_exptime/gain
        hdu[ext].data = (sci-1024.)/spec_exptime/spec_gain
    
    
        ext = 2
        chi2 = hdu[ext].data
        hdu[ext].data = np.sqrt(chi2)/spec_exptime/spec_gain
    
        hdu.writeto(new_slitless, overwrite=True, output_verify='fix')
        print("Writing",new_slitless)
    
    
    # In[ ]:

import datetime
from astropy.time import Time


print(all_slitless)


hdu = pyfits.open(all_slitless[0][0])

#print(hdu[0].header)
#print(hdu[1].header)
#exptime = hdu[0].header['EXPTIME']
crval1 = hdu[1].header['CRVAL1']
crval2 = hdu[1].header['CRVAL2']
#date = hdu[1].header['DATE'] # (YYYY-MM-DDThh:mm:ss UT)
date = '2023-08-09T16:06:51'

time_obj = datetime.datetime.fromisoformat(date)
date_obs = time_obj.strftime("%Y-%m-%d")
time_obs = time_obj.strftime("%H:%M:%S")

t = Time(date, format='fits')
mjd = t.mjd



slitless_dict = {'FILTER':'RED',
                 'INSTRUME':'NISP-GLWv1',
                 'TARGNAME':'Euclid',
                 'DATE-OBS':date_obs ,
                 'TIME-OBS':time_obs, 
                 'EXPSTART':mjd,
                 'EXPTIME':spec_exptime,
                 'RA_TARG':crval1, 
                 'DEC_TARG':crval2, 
                 'POSTARG1':0.000000, 
                 'POSTARG2':0.000000,
                 'DETECTOR':'IR',
                 'EXTVER':1,
                 'ABZP':26.285,
                }


roll_ang = 0.0

###############
### SPECTRA ###
###############
for i,slitless_file in enumerate(all_slitless[0]):
    print(slitless_file)


    
    new_slitless_file = slitless_file.replace(".fits",suffix+".fits")
    all_final_slitless.append(new_slitless_file)
        
    hdu = pyfits.open(slitless_file)
    hdu.info()

    if all_zodi:
        zodi_file = all_zodi[i]
        print(zodi_file)
        hdu_zodi = pyfits.open(zodi_file)
        zodi = hdu_zodi[ext].data
    else:
        zodi = 0.0
#for i,sf in enumerate(all_slitless[0]):    

    
#    hdu = pyfits.open(sf)

    det_id = hdu[1].header['DET_ID']
    
    exptime = slitless_dict['EXPTIME']
    
    #roll_ang = float(sf.split("_")[4])
    slitless_dict['PA_V3'] = roll_ang
    dtime = exptime/(60.*60.*24.)
    slitless_dict['EXPSTART'] += dtime*i

    hdu1 = add_grizli_headers(hdu, slitless_dict)

    ext = 1
    sci = hdu1[ext].data

    #hdu[ext].data = sci/spec_exptime/gain
    #hdu1[ext].data = (sci-1024.)/spec_exptime/spec_gain
    hdu1[ext].data = (1.0*sci - 1.0*zodi)/spec_exptime/spec_gain
    
    
    ext = 2
    chi2 = hdu1[ext].data
    hdu1[ext].data = np.sqrt(chi2)/spec_exptime/spec_gain
    
    #hdu1.writeto(new_slitless, overwrite=True, output_verify='fix')
    #print("Writing",new_slitless)
    
    #nsf = "-".join([prefix,"%03d.0" % float(roll_ang),grism]) + "_flt.fits"
    new_slitless = "-".join([prefix,"det",det_id,"%03d.0" % float(roll_ang),grism]) + "_flt.fits"
    hdu1.writeto(new_slitless, overwrite=True, output_verify='fix')
    print("Writing",new_slitless)


#sys.exit()


nslitless_files = glob.glob('*-red_flt.fits')
print(nslitless_files)




    
    
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

#for slitless in all_final_slitless:
for slitless in nslitless_files:
    wcs_pixel_scale(slitless)


# In[ ]:

def flux2mag(flux, fluxerr=None):
    """
    flux in uJy
    mag in AB
    """
    mag = -2.5*np.log10(flux) + 23.9
    if fluxerr is not None:
        magerr = 2.5*(fluxerr/(flux*np.log(10)))
        return mag, magerr
    else: 
        return mag
        
print(all_cat[0])
phot = Table.read(all_cat[0], format="fits") # ref_cat in multimission
print(all_cat[0]," number_of_sources =",len(phot))

print(phot.colnames)
#print(phot)
#print(phot["vis_flux_aper_0","vis_flux_aper_1","vis_flux_aper_2","vis_flux_aper_3","vis_flux_aper_4","vis_flux_aper_5","vis_flux_aper_6"])
#print(phot["nisp_y_flux_aper_0","nisp_y_flux_aper_1","nisp_y_flux_aper_2","nisp_y_flux_aper_3","nisp_y_flux_aper_4","nisp_y_flux_aper_5","nisp_y_flux_aper_6"])
#print(phot["nisp_j_flux_aper_0","nisp_j_flux_aper_1","nisp_j_flux_aper_2","nisp_j_flux_aper_3","nisp_j_flux_aper_4","nisp_j_flux_aper_5","nisp_j_flux_aper_6"])
#print(phot["nisp_h_flux_aper_0","nisp_h_flux_aper_1","nisp_h_flux_aper_2","nisp_h_flux_aper_3","nisp_h_flux_aper_4","nisp_h_flux_aper_5","nisp_h_flux_aper_6"])


bands = ["vis","nisp_y","nisp_j","nisp_h"]

#aper = range(6+1)

#all_keys = []
#for b in bands:
#    for a in aper:
#        key = "_".join([b,"flux","aper","%i" % (a)])
#        print(key)
#        #all_keys.append()
#        phot[key].info.format = "%.4f"



def curve_of_growth(flux, size, fluxerr=None):

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    if fluxerr is not None:
        ax1.errorbar(size, flux, yerr=fluxerr)
    else:
        ax1.scatter(size, flux)

    plt.show()


import astropy.units as u
from grizli.prep import SEXTRACTOR_PHOT_APERTURES

print("phot apertures (pixels) =", SEXTRACTOR_PHOT_APERTURES)

pixel_scale = 0.06
#pixel_scale = 0.3
SEXTRACTOR_PHOT_APERTURES_ARCSEC = [float(ap)*pixel_scale for ap in SEXTRACTOR_PHOT_APERTURES.split(',')]
#SEXTRACTOR_PHOT_APERTURES_ARCSEC = [float(ap)*pixel_scale*u.arcsec for ap in SEXTRACTOR_PHOT_APERTURES.split(',')]
print("phot apertures (arcsec) =", SEXTRACTOR_PHOT_APERTURES_ARCSEC)



idx = 0
flux = np.array(list(phot[["vis_flux_aper_{}".format(i) for i in range(6+1)]][idx]))
fluxerr = np.array(list(phot[["vis_fluxerr_aper_{}".format(i) for i in range(6+1)]][idx]))
#print(np.array(list(phot[["vis_mask_aper_{}".format(i) for i in range(6+1)]][idx])))
#print(np.array(list(phot[["vis_bkg_aper_{}".format(i) for i in range(6+1)]][idx])))
#print(np.array(list(phot[["vis_flag_aper_{}".format(i) for i in range(6+1)]][idx])))

print(flux)
print(fluxerr)
mag, magerr = flux2mag(flux, fluxerr=fluxerr)
print(mag)
size = SEXTRACTOR_PHOT_APERTURES_ARCSEC
curve_of_growth(flux, size, fluxerr=fluxerr)
curve_of_growth(mag, size, fluxerr=magerr)

#print(phot["vis_flux_aper_0","vis_flux_aper_1","vis_flux_aper_2","vis_flux_aper_3"])
#print(phot["vis_flux_aper_0","vis_flux_aper_1","vis_flux_aper_2","vis_flux_aper_3"].as_array())
#print(flux2mag(np.array(phot["vis_flux_aper_0"])))
#print(flux2mag(np.array(phot["vis_flux_aper_1"])))
#print(flux2mag(np.array(phot["vis_flux_aper_2"])))
#print(flux2mag(np.array(phot["vis_flux_aper_3"])))
#print(flux2mag(np.array(phot["vis_flux_aper_4"])))
#print(flux2mag(np.array(phot["vis_flux_aper_5"])))
#print(flux2mag(np.array(phot["vis_flux_aper_6"])))
#print(phot["nisp_y_flux_aper_0","nisp_y_flux_aper_1","nisp_y_flux_aper_2","nisp_y_flux_aper_3"])
#print(phot["nisp_j_flux_aper_0","nisp_j_flux_aper_1","nisp_j_flux_aper_2","nisp_j_flux_aper_3"])
#print(phot["nisp_h_flux_aper_0","nisp_h_flux_aper_1","nisp_h_flux_aper_2","nisp_h_flux_aper_3"])

print("mag_auto =", phot["mag_auto"][idx])
print("flux_auto =", phot["flux_auto"][idx])
print("kron_rcirc =", phot["kron_rcirc"][idx])
print("kron_radius =", phot["kron_radius"][idx])
print(phot.meta["MINKRON"])


print(bands)


# add magnitudes to photometry
for b in bands:
    for a in range(6+1):
        flux_key = "{0}_flux_aper_{1}".format(b,a)
        fluxerr_key = "{0}_fluxerr_aper_{1}".format(b,a)
        mag_key = "{0}_mag_aper_{1}".format(b,a)
        magerr_key = "{0}_magerr_aper_{1}".format(b,a)

        flux = phot[flux_key]
        fluxerr = phot[fluxerr_key]

        mag, magerr = flux2mag(flux,fluxerr=fluxerr)

        phot[mag_key] = mag
        phot[magerr_key] = magerr
        



bands_aper = ["{0}_mag_aper_{1}".format(b,aper) for b in bands]
bands_flux_aper = ["{0}_flux_aper_{1}".format(b,aper) for b in bands]



# In[ ]:

if plot:
    fig = plt.figure(figsize=(6,6))
    
    k = 0

    for b in bands_aper:


          

        print(np.nanmin(phot[b]),np.nanmax(phot[b]))
        
        ax = fig.add_subplot(2,2,k+1)
        ax.hist(phot[b],range=(10,32),bins=44)
    
    
        ax.set_title(b, fontsize=12)
        ax.set_xlabel("mag")
        ax.set_ylabel("N")
        ax.set_xlim(10,32)
        
        k += 1

        
    plt.show()
    
    



print(all_ref_files)

#for phot in all_phot:
#    print(phot)






if plot:

    
    fig = plt.figure(figsize=(6,6))
    
    k = 0
    #for direct,phot in zip(all_ref_files,all_phot):
    for direct,b in zip(all_ref_files[:-1],bands_aper):
    
    
        hdu = pyfits.open(direct)
        img = hdu[1].data

        print(b)
        filt = phot[b] < plot_mag_cut
        print(len(phot))

        subphot = phot[filt]
        print(len(subphot))
        print()
        
        ax1 = fig.add_subplot(2,2,k+1)
        ax1.imshow(img,origin="lower", vmin=-0.1, vmax=0.1)
        ax1.scatter(subphot[x_image_key], subphot[y_image_key], fc="None", ec="r", s=50, lw=0.5, alpha=1.0)
    
        if plot_labels:
            for j,src in enumerate(subphot[b]):
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


# remove weird source from catalog
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

#print(all_final_direct[i])
print("Total number of sources found =",len(phot))
print("Search primer and all_tbl for RA/DEC matches")

#c_prime = SkyCoord(ra=primer["RA"], dec=primer["DEC"])
#c_prime = SkyCoord(ra=primer["RA_MAG"], dec=primer["DEC_MAG"])
c_prime = SkyCoord(ra=primer["RA_MAG"]*u.deg, dec=primer["DEC_MAG"]*u.deg)
c_phot = SkyCoord(ra=phot[x_world_key], dec=phot[y_world_key])

#idx, d2d, d3d = c_prime.match_to_catalog_sky(c_phot)
idx, d2d, d3d = c_phot.match_to_catalog_sky(c_prime)

filt = d2d < 1*u.arcsec

if plot:
    ax1 = fig.add_subplot(111)
    ax1.hist(d2d[filt].value*3600.,bins=25)
    ax1.set_xlabel("Distance [arcsec]")
    #ax1.set_title(all_ref_files[i], fontsize=12)

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
#print(list(match_clean_tbl['idx']))
print("Number of matches =",len(match_clean_tbl))
print()

if plot:
    plt.show()
        





if plot:
    
    fig = plt.figure(figsize=(10,6))
    
    k = 0    
    for direct,b in zip(all_ref_files[:-1],bands_aper):
    
    
        hdu = pyfits.open(direct)
        img = hdu[1].data
    
        ax1 = fig.add_subplot(2,2,k+1)
        ax1.imshow(img,origin="lower",vmin=-0.1,vmax=0.1)
        #ax1.scatter(phot[x_image_key],phot[y_image_key],fc="None",ec="r",s=60)


        filt = match_clean_tbl[b] < plot_mag_cut
        
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



    
    
#print("N =",len(all_phot_matched_clean))
#print()

#print(all_phot_matched_clean)

verb = 0

for i,b in enumerate(bands_flux_aper):
    zps = []
    mags = []

    
    
    for k,row in enumerate(match_clean_tbl):
        
        #ind = np.argmax(match_clean_tbl[flux_key])

        #print("ind =",ind)
        dn = row[b]
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

#sys.exit()

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
#print(match_clean_tbl)
print(match_clean_tbl[-1][id_key, x_image_key, y_image_key, mag_key])

T_END = time.time()

print()
print("Finished in %.1f seconds" % (T_END-T_START))
#sys.exit()

with open('match_clean_tbl.pickle', 'wb') as f:
    # Pickle the 'data' dictionary using the highest protocol available.
    pickle.dump(match_clean_tbl, f, pickle.HIGHEST_PROTOCOL)

yaml_dict["all_final_slitless"] = all_final_slitless

os.chdir(YAML_PATH)
with open(os.path.basename(yaml_file), 'w',) as f:
    yaml.dump(yaml_dict, f, sort_keys=False)


