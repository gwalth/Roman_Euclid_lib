
import glob

from astropy.table import Table
import astropy.io.fits as pyfits
import astropy.wcs as pywcs

import numpy as np
from scipy import integrate, interpolate

import matplotlib.pyplot as plt


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

#print(cat.colnames)
#print(len(cat))


points = np.array([[a,d] for a,d in cat["X_WORLD","Y_WORLD"]])
#print(points)

ext = 1
slitless_file = glob.glob(simdata_path + "Euclid_DET11_slitless.fits")[0]
#print(slitless_file)
pf = pyfits.open(slitless_file, ignore_missing_end=1)
head = pf[ext].header
#img = pf[ext].data
pw = pywcs.WCS(fobj=pf, header=head)
A = pw.calc_footprint()
#print(A)
#inside_det11 = []
#for p0 in points:
#    #print(p0)
#    inside = insidePolygon(A, p0)
#    #print(inside)
#
#    inside_det11.append(inside)

inside_det11 = np.array([insidePolygon(A, p0) for p0 in points]) 
inside_det11 = np.multiply(inside_det11, 1)

#########################################
# return indices of filter (boolean mask)
#########################################
det11_indices = np.where(inside_det11)[0]
#print(det11_indices)
#cat_ids = [cat['NUMBER'][i] for i in det11_indices]
#print(cat_ids)



#points_culled = points[inside_det11]    
points_culled = points[det11_indices]    
#print(points_culled)
#print(len(points_culled))
#skycoords = pw.wcs_pix2world(pixcoords,1)
pixcoords = pw.wcs_world2pix(points_culled,1)
#print(pixcoords)
#xc, yc = pixcoords[0]



det11_mag_idx = [[i,m] for i,m in enumerate(cat[det11_indices]["MAG_AUTO"])]
det11_mag_idx_sorted = sorted(det11_mag_idx, key=lambda obj: obj[1])
det11_idx = [idx for idx,mag in det11_mag_idx_sorted]
#print(det11_mag_idx_sorted)
#print(len(det11_mag_idx_sorted))


#print(points_culled[det11_idx])
#print(pixcoords[det11_idx])

#Nfirst = 3
Nfirst = 20
#Nfirst = -1 # all

#ratio_spectra = []
lw = 0.8
alpha = 0.6


verb = 1

fig1 = plt.figure(figsize=(12,4))

ax1 = fig1.add_subplot(131)
ax2 = fig1.add_subplot(132)
ax3 = fig1.add_subplot(133)

#fig2 = plt.figure(figsize=(8,8))
#rows = 5
#cols = 4



for p0, c0 in zip(points_culled[det11_idx][:Nfirst], pixcoords[det11_idx][:Nfirst]):
    #euc_dat.find_source(ra, dec)
    #euc_dat.get_spectra_ids()

    a0, d0 = p0
    xc, yc = c0

    model_wav, model_flux = euc_dat.get_model_spectra(a0, d0)
    #euc_dat.plot_model_spectra()


    aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=0.0)
    #aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=-0.44)

    aXe.extract_cutout()

    if verb > 1:
        aXe.plot_cutout(vmin=-0.075, vmax=0.075, verb=1)
        #aXe.plot_trace()
        #aXe.plot_trace(x0=1350, x1=2350, y0=2800, y1=3800)
    aXe.get_thumbnail(a0, d0)

    if verb > 1:
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
        
        #aXe.plot_footprint()
        
        #aXe.plot_sources_wcs(
        #    x=subcat["X_WORLD"], 
        #    y=subcat["Y_WORLD"], 
        #    labels=subcat["MAG_AUTO"], 
        #    vmin=-1, 
        #    vmax=500,
        #    plot_all=True,
        #    mag_cut=17,
        #)
        
        #aXe.plot_extraction()
        
        #y0 = -1e-16
        #y1 = 1e-15
        #y1 = 1e-13
        
        #aXe.plot_extraction(flux_units=True,)# y0=y0, y1=y1)
        
        #aXe.plot_extraction(flux_units=True, y0=y0, y1=y1, #log=True,  
        #                    model=[model_wav, model_flux])
  
    if verb > 1:
        aXe.plot_extraction(flux_units=True,
                            model=[model_wav, model_flux])
        
        #help(EuclidFilters)
        ef = EuclidFilters()
        #ef.plot_filters()
        #print(ef.filters.keys())
        model_fnu = flam_to_fnu(flam=model_flux, wav=model_wav)
        ef.calc_spectra_mag("NISP_H", model_fnu, model_wav/10.)

    interp_model = interpolate.interp1d(model_wav, model_flux)

    flux = aXe.spectrum[2]
    err = aXe.spectrum[3]
    wav = aXe.wav



    ax1.plot(model_wav,model_flux, lw=lw, alpha=alpha)

    ax2.plot(wav, flux, lw=lw, alpha=alpha)

    flux_ratio = np.array([interp_model(w0) for w0 in wav])/flux

    ax3.plot(wav, flux_ratio, lw=lw, alpha=alpha)




    #ratio_spectra.appenda([flux_ratio, wav])



ax1.set_xlim(13600,18700)
ax1.set_ylim(3e-18, 6e-16)
ax1.set_yscale("log")
ax1.set_title("Input Galaxy Models")

ax2.set_xlim(13600,18700)
ax2.set_ylim(3e-18, 6e-16)
ax2.set_yscale("log")
ax2.set_title("Extracted Spectra")

ax3.set_xlim(13600,18700)
ax3.set_ylim(-0.5,4)
ax3.set_title("Input Galaxy Model/Extracted spectra [ratio]")

plt.show()
