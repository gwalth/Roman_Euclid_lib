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

aper = 20.
#aper = 8.
#aper = 16.

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
#aper = 8.

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







points = np.array([[a,d] for a,d in cat1["RA_MAG","DEC_MAG"]])


print(points)
print(len(points))


ext = 1
slitless_file = glob.glob(simdata_path + "Euclid_DET11_slitless.fits")[0]
#print(slitless_file)
pf = pyfits.open(slitless_file, ignore_missing_end=1)
head = pf[ext].header
#img = pf[ext].data
pw = pywcs.WCS(fobj=pf, header=head)
A = pw.calc_footprint()


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



det11_mag_idx = [[i,m] for i,m in enumerate(cat1[det11_indices]["NIR_H"])]
det11_mag_idx_sorted = sorted(det11_mag_idx, key=lambda obj: obj[1])
det11_idx = [idx for idx,mag in det11_mag_idx_sorted]
#print(det11_mag_idx_sorted)
#print(len(det11_mag_idx_sorted))


#print(points_culled[det11_idx])
#print(pixcoords[det11_idx])

#Nfirst = 3
#Nfirst = 20
Nfirst = 30
#Nfirst = -1 # all

#N0 = 50
#N1 = 60

#ratio_spectra = []
lw = 0.8
alpha = 0.6


verb = 0

fig1 = plt.figure(figsize=(12,4))

ax1 = fig1.add_subplot(131)
ax2 = fig1.add_subplot(132)
ax3 = fig1.add_subplot(133)

#fig2 = plt.figure(figsize=(8,8))
#rows = 5
#cols = 4

for p0, c0 in zip(points_culled[det11_idx][:Nfirst], pixcoords[det11_idx][:Nfirst]):
#for p0, c0 in zip(points_culled[det11_idx][N0:N1], pixcoords[det11_idx][N0:N1]):

    a0, d0 = p0
    xc, yc = c0

    dr = np.sqrt((cat1['RA_MAG']-a0)**2*np.cos(d0/180*np.pi)**2 + 
         (cat1['DEC_MAG']-d0)**2)*3600.
    id = cat1['SOURCE_ID'][np.argmin(dr)]
    ind = np.argmin(dr)
    obj_mag = cat1['NIR_H'][np.argmin(dr)]
    print('ID:%d, IND:%d, mag=%.2f, dr=%.2f"' % (id, ind, obj_mag, np.min(dr)))

    model_wav, model_flux = get_model_spectra_A(model_spectra, ind)
    if verb:
        plot_model_spectra_A(model_spectra, ind)



    aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=0.0)
    #aXe.get_beam_cutout(xc, yc, aper, beam='A', offset=-0.44)

    aXe.extract_cutout()

    if verb > 1:
        aXe.plot_cutout(vmin=-0.075, vmax=0.075, verb=1)
        #aXe.plot_trace()
        #aXe.plot_trace(x0=1350, x1=2350, y0=2800, y1=3800)
    aXe.get_thumbnail(a0, d0)

    if verb > 1:
        aXe.plot_thumbnail(cat=cat1, ra_key="RA_MAG", dec_key="DEC_MAG", 
                       label_key="NIR_H")
        
        #aXe.plot_footprint()
        
        #aXe.plot_extraction()
        
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




ax1.set_xlim(13600,18700)
#ax1.set_ylim(3e-18, 6e-16)
ax1.set_yscale("log")
ax1.set_title("Input Galaxy Models")

ax2.set_xlim(13600,18700)
#ax2.set_ylim(3e-18, 6e-16)
ax2.set_yscale("log")
ax2.set_title("Extracted Spectra")

ax3.set_xlim(13600,18700)
ax3.set_ylim(-0.5,4)
ax3.set_title("Input Galaxy Model/Extracted spectra [ratio]")

plt.show()
