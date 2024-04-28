import glob

from astropy.table import Table
import astropy.io.fits as pyfits
import astropy.wcs as pywcs

import numpy as np
from scipy import integrate, interpolate

import matplotlib.pyplot as plt


#import grizli_aXe_glw
from grizli_aXe_glw import aXeConf_glw, EuclidData, EuclidFilters, zodi_subtract, insidePolygon, flam_to_fnu, get_model_spectra_A, plot_model_spectra_A


#simdata_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/TestPoints/Prep/"
simdata_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/SIM_10_18_22/Prep/"
raw_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/Raw/EuclidSIMS/testARRAY_2023_06_23/frame1/"

sci_file = raw_path  + "NISPS_TIPS_testARRAY_2023_06_23_frame1.fits"
zodi_file = raw_path  + "NISPS_TIPS_testARRAY_JustZodi_2023_06_23_frame1.fits"
out_file = simdata_path  + "NISPS_TIPS_testARRAY_subZodi_2023_06_23_frame1.fits"


# In[ ]:


zodi_subtract(sci_file, zodi_file, out_file, spec_gain=6.0)
#zodi_subtract(sci_file, zodi_file, out_file, spec_gain=1.0)
#zodi_subtract(sci_file, zodi_file, out_file, spec_gain=2.0)


# In[ ]:


#prefix = "/Users/gwalth/data/aXeSIM/Roman/aXeSIM_Roman/"

conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/"        # BEAMA
#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/CONF11/"  # BEAM_A
#roman_path = prefix + "roman_setup/"
simdata_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/SIM_10_18_22/Prep/"

conf = "NISP_RGS000_11.conf"
sens_file = 'SENS_A.fits'
#input_spec_list = "input_spectra.lis"

direct_file = simdata_path + "Euclid_NISP_H_ref.fits"
#slitless_file = simdata_path + "GRS_FOV3_roll0_dx0_dy0_SCA1_slitless_final.fits"
#slitless_file = simdata_path + "Euclid_DET11_slitless.fits"
#slitless_file = simdata_path + "Euclid-det-11-000.0-red_flt.fits"
slitless_file = simdata_path + "NISPS_TIPS_testARRAY_subZodi_2023_06_23_frame1.fits"


aper = 8.

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
aXe.load_files(direct_file, slitless_file, offset=0)
#aXe.load_files(direct_file, slitless_file, offset=-1024.)

aXe.load_sensitivity(conf_path+sens_file)
aper = 8.

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
aXe.load_files(direct_file, slitless_file, offset=0)
#aXe.load_files(direct_file, slitless_file, offset=-1024.)

aXe.load_sensitivity(conf_path+sens_file)


# In[ ]:


catalog1 = "CATALOG_WP9_INTERMEDIATE_RUN_v2_NewCoord_mod.fits"
catalog2 = "EUC_SIM_NIS-TEST_1717_catalog_11.cat"

catalog1_file = simdata_path + catalog1
catalog2_file = simdata_path + "../../Raw/SIM_10_18_22_singleframe/data/" + catalog2

##########
# catalog1
##########
cat1 = Table.read(catalog1_file, format="fits")
N1 = len(cat1)
cat1["index"] = np.arange(N1)
print([col for col in cat1.colnames if "TU_" in col])
        
Euclid_bands = ['VIS','NIR_Y','NIR_J','NIR_H']
Euclid_bands_flux = ['TU_FNU_VIS_MAG', 'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 
                                  'TU_FNU_H_NISP_MAG'] 

for bk,fk in zip(Euclid_bands, Euclid_bands_flux):
    fnu_Jy = cat1[fk] # Jy    
    mab = -2.5*np.log10(fnu_Jy) + 8.90   # Jy --> AB mag            
    mab[np.isinf(mab)]=-99.
    cat1[bk] = mab  
    
print(cat1.colnames)
print(len(cat1))

##########
# catalog2
##########
cat2 = Table.read(catalog2_file, format="ascii") 
N2 = len(cat2)
cat2["index"] = np.arange(N2)
print(N2)


# In[ ]:


model_spectra = simdata_path + "A.fits"
print(model_spectra)


# ## Source 1

# In[ ]:


########################################################
# find a source in all of the Euclid detector footprints
########################################################
# bright source near the center of DET11
a0, d0 = 228.6881592, 6.3239199

dr = np.sqrt((cat1['RA_MAG']-a0)**2*np.cos(d0/180*np.pi)**2 + 
             (cat1['DEC_MAG']-d0)**2)*3600.
id = cat1['SOURCE_ID'][np.argmin(dr)]
ind = np.argmin(dr)
obj_mag = cat1['NIR_H'][np.argmin(dr)]
print('ID:%d, IND:%d, mag=%.2f, dr=%.2f"' % (id, ind, obj_mag, np.min(dr)))

#############################################################
# conversion between the reference image and slitless spectra
#############################################################

import glob

import astropy.wcs as pywcs

slitless_files = glob.glob(simdata_path + "Euclid_DET??_slitless.fits")
#slitless_files = [simdata_path + "Euclid_DET11_slitless.fits"]
slitless_files.sort()
print(slitless_files)

# find a single coordinate in all slitless DET images
for sf in slitless_files:

    ext = 1
    print(sf)
    pf = pyfits.open(sf, ignore_missing_end=1)
    head = pf[ext].header
    #img = pf[ext].data
    pw = pywcs.WCS(fobj=pf, header=head)
    A = pw.calc_footprint()
    print(A)

    p0 = [a0, d0]
    print(p0)
    inside = insidePolygon(A, p0)
    print(inside)
    
    if inside:
        #skycoords = pw.wcs_pix2world(pixcoords,1)
        pixcoords = pw.wcs_world2pix([p0],1)
        print(pixcoords)
        xc, yc = pixcoords[0]


# In[ ]:


#euc_dat.find_source(ra, dec)
#euc_dat.get_spectra_ids()
model_wav, model_flux = get_model_spectra_A(model_spectra, ind)
plot_model_spectra_A(model_spectra, ind)


# In[ ]:


aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=0.0)
#aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=-0.38)

aXe.extract_cutout()

aXe.plot_cutout(vmin=-0.5, vmax=0.5, verb=1)
#aXe.plot_trace()
#aXe.plot_trace(x0=1350, x1=2350, y0=2800, y1=3800)

#aXe.plot_sources(primer, vmin=-1, vmax=1)

aXe.plot_footprint()

#aXe.plot_sources_wcs(
#    x=subcat["X_WORLD"], 
#    y=subcat["Y_WORLD"], 
#    labels=subcat["MAG_AUTO"], 
#    vmin=-1, 
#    vmax=500,
#    plot_all=True,
#    mag_cut=17,
#)

aXe.plot_extraction()

y0 = -1e-15
y1 = 3e-15

#y1 = 1e-13

aXe.plot_extraction(flux_units=True,)# y0=y0, y1=y1)

aXe.plot_extraction(flux_units=True, y0=y0, y1=y1, #log=True,  
                    model=[model_wav, model_flux])


# ## Source 2

# In[ ]:


########################################################
# find a source in all of the Euclid detector footprints
########################################################
# bright source near the bottom of DET11
a0, d0 = 228.6881780, 6.2909390

dr = np.sqrt((cat1['RA_MAG']-a0)**2*np.cos(d0/180*np.pi)**2 + 
             (cat1['DEC_MAG']-d0)**2)*3600.
id = cat1['SOURCE_ID'][np.argmin(dr)]
ind = np.argmin(dr)
obj_mag = cat1['NIR_H'][np.argmin(dr)]
print('ID:%d, IND:%d, mag=%.2f, dr=%.2f"' % (id, ind, obj_mag, np.min(dr)))

#############################################################
# conversion between the reference image and slitless spectra
#############################################################

import glob

import astropy.wcs as pywcs

slitless_files = glob.glob(simdata_path + "Euclid_DET??_slitless.fits")
#slitless_files = [simdata_path + "Euclid_DET11_slitless.fits"]
slitless_files.sort()
print(slitless_files)

# find a single coordinate in all slitless DET images
for sf in slitless_files:

    ext = 1
    print(sf)
    pf = pyfits.open(sf, ignore_missing_end=1)
    head = pf[ext].header
    #img = pf[ext].data
    pw = pywcs.WCS(fobj=pf, header=head)
    A = pw.calc_footprint()
    print(A)

    p0 = [a0, d0]
    print(p0)
    inside = insidePolygon(A, p0)
    print(inside)
    
    if inside:
        #skycoords = pw.wcs_pix2world(pixcoords,1)
        pixcoords = pw.wcs_world2pix([p0],1)
        print(pixcoords)
        xc, yc = pixcoords[0]


# In[ ]:


#euc_dat.find_source(ra, dec)
#euc_dat.get_spectra_ids()
model_wav, model_flux = get_model_spectra_A(model_spectra, ind)
plot_model_spectra_A(model_spectra, ind)


# In[ ]:


aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=0.0)
#aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=-0.38)

aXe.extract_cutout()

aXe.plot_cutout(vmin=-0.5, vmax=0.5, verb=1)
#aXe.plot_trace()
#aXe.plot_trace(x0=1350, x1=2350, y0=2800, y1=3800)

#aXe.plot_sources(primer, vmin=-1, vmax=1)

aXe.plot_footprint()

#aXe.plot_sources_wcs(
#    x=subcat["X_WORLD"], 
#    y=subcat["Y_WORLD"], 
#    labels=subcat["MAG_AUTO"], 
#    vmin=-1, 
#    vmax=500,
#    plot_all=True,
#    mag_cut=17,
#)

aXe.plot_extraction()

y0 = -1e-16
y1 = 3e-16

#y1 = 1e-13

aXe.plot_extraction(flux_units=True,)# y0=y0, y1=y1)

aXe.plot_extraction(flux_units=True, y0=y0, y1=y1, #log=True,  
                    model=[model_wav, model_flux])


# ### Source 3

# In[ ]:


########################################################
# find a source in all of the Euclid detector footprints
########################################################
# bright source near the top of DET11
a0, d0 = 228.6882612, 6.3599378

dr = np.sqrt((cat1['RA_MAG']-a0)**2*np.cos(d0/180*np.pi)**2 + 
             (cat1['DEC_MAG']-d0)**2)*3600.
id = cat1['SOURCE_ID'][np.argmin(dr)]
ind = np.argmin(dr)
obj_mag = cat1['NIR_H'][np.argmin(dr)]
print('ID:%d, IND:%d, mag=%.2f, dr=%.2f"' % (id, ind, obj_mag, np.min(dr)))

#############################################################
# conversion between the reference image and slitless spectra
#############################################################

import glob

import astropy.wcs as pywcs

slitless_files = glob.glob(simdata_path + "Euclid_DET??_slitless.fits")
#slitless_files = [simdata_path + "Euclid_DET11_slitless.fits"]
slitless_files.sort()
print(slitless_files)

# find a single coordinate in all slitless DET images
for sf in slitless_files:

    ext = 1
    print(sf)
    pf = pyfits.open(sf, ignore_missing_end=1)
    head = pf[ext].header
    #img = pf[ext].data
    pw = pywcs.WCS(fobj=pf, header=head)
    A = pw.calc_footprint()
    print(A)

    p0 = [a0, d0]
    print(p0)
    inside = insidePolygon(A, p0)
    print(inside)
    
    if inside:
        #skycoords = pw.wcs_pix2world(pixcoords,1)
        pixcoords = pw.wcs_world2pix([p0],1)
        print(pixcoords)
        xc, yc = pixcoords[0]


# In[ ]:


#euc_dat.find_source(ra, dec)
#euc_dat.get_spectra_ids()
model_wav, model_flux = get_model_spectra_A(model_spectra, ind)
plot_model_spectra_A(model_spectra, ind)


# In[ ]:


aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=0.0)
#aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=-0.38)

aXe.extract_cutout()

aXe.plot_cutout(vmin=-0.5, vmax=0.5, verb=1)
#aXe.plot_trace()
#aXe.plot_trace(x0=1350, x1=2350, y0=2800, y1=3800)

#aXe.plot_sources(primer, vmin=-1, vmax=1)

aXe.plot_footprint()

#aXe.plot_sources_wcs(
#    x=subcat["X_WORLD"], 
#    y=subcat["Y_WORLD"], 
#    labels=subcat["MAG_AUTO"], 
#    vmin=-1, 
#    vmax=500,
#    plot_all=True,
#    mag_cut=17,
#)

aXe.plot_extraction()

y0 = -1e-17
y1 = 6e-17

#y1 = 1e-13

aXe.plot_extraction(flux_units=True,)# y0=y0, y1=y1)

aXe.plot_extraction(flux_units=True, y0=y0, y1=y1, #log=True,  
                    model=[model_wav, model_flux])


