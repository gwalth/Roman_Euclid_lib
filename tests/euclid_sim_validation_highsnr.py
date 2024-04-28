
from astropy.table import Table
import astropy.io.fits as pyfits
import numpy as np

#import grizli_aXe_glw
from grizli_aXe_glw import aXeConf_glw, EuclidData, EuclidFilters, zodi_subtract, insidePolygon, flam_to_fnu

# ## highSNR
# [top](#Table-of-Contents)



simdata_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/TestPoints/Prep/"
raw_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/Raw/EuclidSIMS/TestPoints_highSNR_mod1_14324_2023_05_26/frame1/"

sci_file = raw_path  + "NISPS_TIPS_TestPoints_highSNR_mod1_14324_2023_05_26_frame1.fits"
zodi_file = raw_path  + "NISPS_TIPS_TestPoints_highSNR_JustZodi_14324_2023_05_26_frame1.fits"
out_file = simdata_path  + "NISPS_TIPS_TestPoints_highSNR_subZodi_14324_2023_05_26_frame1.fits"


# In[14]:


zodi_subtract(sci_file, zodi_file, out_file, spec_gain=6.0)
#zodi_subtract(sci_file, zodi_file, out_file, spec_gain=1.0)
#zodi_subtract(sci_file, zodi_file, out_file, spec_gain=2.0)


# In[86]:


#prefix = "/Users/gwalth/data/aXeSIM/Roman/aXeSIM_Roman/"

conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/"        # BEAMA
#conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/CONF11/"  # BEAM_A
#roman_path = prefix + "roman_setup/"
simdata_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/TestPoints/Prep/"

conf = "NISP_RGS000_11.conf"
sens_file = 'SENS_A.fits'
#input_spec_list = "input_spectra.lis"

direct_file = simdata_path + "Euclid-NISP_H_ref.fits"
#slitless_file = simdata_path + "GRS_FOV3_roll0_dx0_dy0_SCA1_slitless_final.fits"
#slitless_file = simdata_path + "Euclid_DET11_slitless.fits"
#slitless_file = simdata_path + "Euclid-det-11-000.0-red_flt.fits"
slitless_file = simdata_path + "NISPS_TIPS_TestPoints_highSNR_subZodi_14324_2023_05_26_frame1.fits"

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


# In[47]:


euc_dat = EuclidData(
        catalog1 = "TestPoints_highSNR_mod1_14324.fits",
        catalog2 = "NIS_detector_catalog_11.cat",
        model_spectra = "NIS_catalog_file_11.spc.fits",
        euclud_sims_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/Raw/EuclidSIMS/",
        catalog1_path = "./",
        catalog2_path = "TestPoints_highSNR_mod1_14324_2023_05_26/frame1/Catalogs/",
        model_spectra_path = "TestPoints_highSNR_mod1_14324_2023_05_26/frame1/Spectra/",
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

print(cat.colnames)
print(len(cat))


# ### Source 1

# In[49]:


########################################################
# find a source in all of the Euclid detector footprints
########################################################
# bright source near the center of DET11
a0, d0 = 214.6623540, 57.8587647

dr = np.sqrt((cat['X_WORLD']-a0)**2*np.cos(d0/180*np.pi)**2 + 
             (cat['Y_WORLD']-d0)**2)*3600.
id = cat['NUMBER'][np.argmin(dr)]
print(id)
print(np.argmin(dr))
obj_mag = cat['MAG_AUTO'][np.argmin(dr)]
print('ID:%d, mag=%.2f, dr=%.2f"' % (id, obj_mag, np.min(dr)))

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


# In[50]:


#euc_dat.find_source(ra, dec)
#euc_dat.get_spectra_ids()
model_wav, model_flux = euc_dat.get_model_spectra(a0, d0)
euc_dat.plot_model_spectra()


# In[87]:


aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=0.0)
#aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=-0.44)

aXe.extract_cutout()

aXe.plot_cutout(vmin=-0.5, vmax=0.5, verb=1)
#aXe.plot_trace()
#aXe.plot_trace(x0=1350, x1=2350, y0=2800, y1=3800)
aXe.get_thumbnail(a0, d0)

aXe.plot_thumbnail(cat=cat, vmin=-1, vmax=1)
#aXe.plot_sources(vmin=-1, vmax=1)
#aXe.plot_sources(
#    x=cat["X_WORLD"], 
#    y=cat["Y_WORLD"], 
#    labels=cat["MAG_AUTO"], 
#    vmin=-1, 
#    vmax=1,
#    plot_all=True,
#    mag_cut=17,
#)

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
y1 = 1e-15

#y1 = 1e-13

aXe.plot_extraction(flux_units=True,)# y0=y0, y1=y1)

aXe.plot_extraction(flux_units=True, y0=y0, y1=y1, #log=True,  
                    model=[model_wav, model_flux])


# In[23]:


#help(EuclidFilters)
ef = EuclidFilters()
ef.plot_filters()
print(ef.filters.keys())
model_fnu = flam_to_fnu(flam=model_flux, wav=model_wav)
ef.calc_spectra_mag("NISP_H", model_fnu, model_wav/10.)


# In[25]:


#inspect_direct(a0, d0, cat, direct_file, ra_key="X_WORLD", dec_key="Y_WORLD", mag_key="MAG_AUTO")


# ### Source 2

# In[ ]:


########################################################
# find a source in all of the Euclid detector footprints
########################################################
# another bright source near the center of DET11
a0, d0 = 214.6360185, 57.9337266

dr = np.sqrt((cat['X_WORLD']-a0)**2*np.cos(d0/180*np.pi)**2 + 
             (cat['Y_WORLD']-d0)**2)*3600.
id = cat['NUMBER'][np.argmin(dr)]
obj_mag = cat['MAG_AUTO'][np.argmin(dr)]
print('ID:%d, mag=%.2f, dr=%.2f"' % (id, obj_mag, np.min(dr)))

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
model_wav, model_flux = euc_dat.get_model_spectra(a0, d0)
euc_dat.plot_model_spectra()


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
y1 = 2e-16

#y1 = 1e-13

aXe.plot_extraction(flux_units=True,)# y0=y0, y1=y1)

aXe.plot_extraction(flux_units=True, y0=y0, y1=y1, #log=True,  
                    model=[model_wav, model_flux])


# ### Source 3

# In[88]:


########################################################
# find a source in all of the Euclid detector footprints
########################################################
# another bright source near the bottom
a0, d0 = 214.7283007, 57.8219983

dr = np.sqrt((cat['X_WORLD']-a0)**2*np.cos(d0/180*np.pi)**2 + 
             (cat['Y_WORLD']-d0)**2)*3600.
id = cat['NUMBER'][np.argmin(dr)]
obj_mag = cat['MAG_AUTO'][np.argmin(dr)]
print('ID:%d, mag=%.2f, dr=%.2f"' % (id, obj_mag, np.min(dr)))

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


# In[89]:


#euc_dat.find_source(ra, dec)
#euc_dat.get_spectra_ids()
model_wav, model_flux = euc_dat.get_model_spectra(a0, d0)
euc_dat.plot_model_spectra()


# In[90]:


aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=0.0)
#aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=-0.44)

aXe.extract_cutout()

aXe.plot_cutout(vmin=-0.5, vmax=0.5, verb=1)
#aXe.plot_trace()
#aXe.plot_trace(x0=1350, x1=2350, y0=2800, y1=3800)
aXe.get_thumbnail(a0, d0)

aXe.plot_thumbnail(cat=cat, vmin=-1, vmax=1)
#aXe.plot_sources(vmin=-1, vmax=1)
#aXe.plot_sources(
#    x=cat["X_WORLD"], 
#    y=cat["Y_WORLD"], 
#    labels=cat["MAG_AUTO"], 
#    vmin=-1, 
#    vmax=1,
#    plot_all=True,
#    mag_cut=17,
#)

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
y1 = 1e-15

#y1 = 1e-13

aXe.plot_extraction(flux_units=True,)# y0=y0, y1=y1)

aXe.plot_extraction(flux_units=True, y0=y0, y1=y1, #log=True,  
                    model=[model_wav, model_flux])


# In[91]:


help(EuclidFilters)
ef = EuclidFilters()
ef.plot_filters()
print(ef.filters.keys())
model_fnu = flam_to_fnu(flam=model_flux, wav=model_wav)
ef.calc_spectra_mag("NISP_H", model_fnu, model_wav/10.)


