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
from grizli_aXe_glw import aXeConf_glw, EuclidData, EuclidFilters, zodi_subtract, insidePolygon, flam_to_fnu
from grizli_functions import euclid_wcs, remove_basic_wcs, stats

import matplotlib as mpl    
mpl.rcParams['font.size'] = 8

# ## highSNR
# [top](#Table-of-Contents)

def gauss_fn(param, x):
    # From Dan's gaussfit
    """Parameters: a0 = height of exp, a1 = center of exp, a2 = sigma
    (the width)"""
    a0, a1, a2, a3 = param
    return a0 * np.exp(-(((x - a1) / a2) ** 2.0) / 2.0) + a3


noise = 1.
err_fn = lambda p, x, z: (gauss_fn(p, x) - z) / noise



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
Nfirst = 18
#Nfirst = 24
#Nfirst = 36
#Nfirst = -1 # all


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



#sci_file = raw_path  + "NISPS_TIPS_TestPoints_highSNR_mod1_14324_2023_05_26_frame1.fits"
#zodi_file = raw_path  + "NISPS_TIPS_TestPoints_highSNR_JustZodi_14324_2023_05_26_frame1.fits"
#out_file = simdata_path  + "NISPS_TIPS_TestPoints_highSNR_subZodi_14324_2023_05_26_frame1.fits"



#zodi_subtract(sci_file, zodi_file, out_file, spec_gain=6.0)
#zodi_subtract(sci_file, zodi_file, out_file, spec_gain=1.0)
#zodi_subtract(sci_file, zodi_file, out_file, spec_gain=2.0)


# In[86]:


#prefix = "/Users/gwalth/data/aXeSIM/Roman/aXeSIM_Roman/"

conf_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv1/"
pre_process_path = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS000p0/"

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
conf = "test3_NISP_RGS000_12.conf"

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



det_mag_idx = [[i,m] for i,m in enumerate(cat[det_indices]["MAG_AUTO"])]
det_mag_idx_sorted = sorted(det_mag_idx, key=lambda obj: obj[1])
det_idx = [idx for idx,mag in det_mag_idx_sorted]
det_mag = [mag for idx,mag in det_mag_idx_sorted]

print(det_mag_idx_sorted)
print(det_idx)
print(det_mag)
print(len(det_mag_idx_sorted))




print(points_culled[det_idx])
print(pixcoords[det_idx])
print(cat[det_indices]["MAG_AUTO"][det_idx])





#ratio_spectra = []

fig1, axes1 = plt.subplots(nrows=1, ncols=3, figsize=(12,4))
ax11, ax12, ax13 = axes1


ncols = 2
nrows = math.ceil(Nfirst / ncols)
fig2, axes2 = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12,6))
fig2.subplots_adjust(wspace=0.0, hspace=0.0)

ncols = 6
nrows = math.ceil(Nfirst / ncols)
fig3, axes3 = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6,6))
fig3.subplots_adjust(wspace=0.0, hspace=0.0)

ncols = 6
nrows = math.ceil(Nfirst / ncols)
fig4, axes4 = plt.subplots(nrows=nrows, ncols=ncols)
fig4.subplots_adjust(wspace=0.0, hspace=0.0)




#plt.subplots_adjust(left=0.1,
#                    bottom=0.1, 
#                    right=0.9, 
#                    top=0.9, 
#                    wspace=0.4, 
#                    hspace=0.4)


#subcat = cat["MAG_AUTO"][det_idx]
#print(subcat)
#print("I am here bitches!")



vmin = np.min(det_mag[:Nfirst])
vmax = np.max(det_mag[:Nfirst])
colors_mag = [cmap( m/(vmax-vmin) - vmin/(vmax-vmin) ) for m in det_mag[:Nfirst]]



radius = cat["KRON_RADIUS"][det_indices][det_idx]

#vmin = np.min(radius[:Nfirst])
#vmax = np.max(radius[:Nfirst])
#colors_rad = [cmap( r/(vmax-vmin) - vmin/(vmax-vmin) ) for r in radius[:Nfirst]]


colors = colors_mag
#colors = colors_rad

#print(det_mag[:Nfirst])
print(vmin)
print(vmax)

normalizer = Normalize(vmin,vmax)
im = cm.ScalarMappable(norm=normalizer, cmap=cmap)


#sys.exit()

sub_radius = radius[:Nfirst]



for i,(p0, c0) in enumerate(zip(points_culled[det_idx][:Nfirst], pixcoords[det_idx][:Nfirst])):


    #euc_dat.get_spectra_ids()

    a0, d0 = p0
    xc, yc = c0

    sid, ind = euc_dat.find_source(a0, d0)
    model_wav, model_flux = euc_dat.get_model_spectra(a0, d0)
    #euc_dat.plot_model_spectra()

    mag = det_mag[i]

    y0 = np.min(model_flux)
    y1 = np.max(model_flux)

    if y0 < np.abs(ymin):
        ymin = y0

    if y1 > ymax:
        ymax = y1



    #aXe.get_beam_invert_cutout(xc, yc, aper, beam='A', offset=0.0)
    aXe.get_beam_invert_cutout(xc, yc, aper, beam='A', offset=0.0, dy0=-300.0, dy1=300.0)
    cutout = aXe.cutout

    ax2 = axes2.flat[i]
    ax2.imshow(cutout[0].T, origin="lower")
    #aXe.get_beam_invert_cutout(xc, yc, aper, beam='A', offset=-0.44)
    ax2.text(0.01, 0.65, "%s" % ind, transform=ax2.transAxes, color="w")
    ax2.text(0.01, 0.10, "%.2f" % mag, transform=ax2.transAxes, color="w")
    if ind in skip_indicies:
        ax2.text(0.5, 0.5, "X", color="r", transform=ax2.transAxes, fontsize=20, ha="center", va="center")



    if verb > 1:
        aXe.plot_cutout(vmin=-0.075, vmax=0.075, verb=1)
        #aXe.plot_trace()
        #aXe.plot_trace(x0=1350, x1=2350, y0=2800, y1=3800)
    aXe.get_thumbnail(a0, d0)
    thumbnail, _ = aXe.thumbnail

    stat_dict = stats(thumbnail.flatten())
    x_std = stat_dict["std"]

    ax3 = axes3.flat[i]
    ax3.imshow(thumbnail, origin="lower", vmin=-0.5*x_std, vmax=5*x_std)
    ax3.axes.get_xaxis().set_visible(False)
    ax3.axes.get_yaxis().set_visible(False)
    ax3.text(0.01, 0.80, "%i" % ind, transform=ax3.transAxes, color="w")
    ax3.text(0.01, 0.10, "%.2f" % mag, transform=ax3.transAxes, color="w")
    if ind in skip_indicies:
        ax3.text(0.5, 0.5, "X", color="r", transform=ax3.transAxes, fontsize=50, ha="center", va="center")


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
  
    #aXe.extract_cutout(disp_axis=1, spatial_axis=0)
    aXe.extract_cutout(disp_axis=1, spatial_axis=0, flip_axis=1)
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
    profile = aXe.profile

    yp_min= np.min(profile)
    yp_max= np.max(profile)
    x_cen = len(profile)/2.
    x_data = np.arange(len(profile))


    ax4 = axes4.flat[i]
    ax4.plot(profile)
    ax4.plot([x_cen,x_cen], [yp_min,yp_max], "--", c="k", alpha=0.5)
    ax4.ticklabel_format(axis='y', style='sci', scilimits=(0, 0)) 
    ax4.text(0.01, 0.80, "%s" % ind, transform=ax4.transAxes)

    if ind in skip_indicies:
        ax4.text(0.5, 0.5, "X", color="r", transform=ax4.transAxes, fontsize=50, ha="center", va="center")

    

    if gfit:

        param = [
                yp_max,
                x_cen,
                len(profile)/3.,
                0,
        ]
        
        result = least_squares(err_fn, param, args=(x_data, profile))
        #print(result)
        
        U, s, Vh = linalg.svd(result.jac, full_matrices=False)
        tol = np.finfo(float).eps * s[0] * max(result.jac.shape)
        w = s > tol
        cov = (Vh[w].T / s[w] ** 2) @ Vh[w]  # robust covariance matrix
        perr = np.sqrt(np.diag(cov))  # 1sigma uncertainty on fitted parameters
        
        dof = result.fun.size - result.x.size
        chi2 = np.sum(result.fun**2)
        rchi2 = chi2 / dof
        
        cov *= rchi2
        perr = np.sqrt(
            np.diag(cov)
        )  # 1sigma uncertainty on fitted parameters
        # print(perr)
        
        print("chisq = %.2f" % chi2)
        # dof is degrees of freedom
        print("degrees of freedom = %.0f" % dof)
        print("reduced chisq = %.2f" % rchi2)
        
        p = result.x
        
        print(p)
        print(perr)
        print()

        x_inc = aper / 100.0
        x_small = np.arange(0, aper, x_inc)
        g1 = gauss_fn(p, x_small)
        ax4.plot(x_small, g1, c="g")






    xl = flux.shape[0]
    #gain = xl*aper


    flux /= gain
    err /= gain

    ax11.plot(model_wav, model_flux, lw=lw, alpha=alpha, color=colors[i])

    ax12.plot(wav, flux, lw=lw, alpha=alpha, color=colors[i])

    flux_ratio = np.array([interp_model(w0) for w0 in wav])/flux

    if ind in skip_indicies:
        ax13.plot(wav, flux_ratio, lw=lw, color="k", alpha=0.2)
    else:
        ax13.plot(wav, flux_ratio, lw=lw, alpha=alpha, color=colors[i])

    print("KRON RADIUS: %.2f" % (sub_radius[i]))
    print()




    #ratio_spectra.appenda([flux_ratio, wav])

print(cat)
print(cat.colnames)


ax11.set_xlim(wmin, wmax)
#ax11.set_ylim(3e-18, 6e-16)
ax11.set_ylim(ymin/2, ymax*2)
ax11.set_yscale("log")
ax11.set_title("Input Galaxy Models")

ax12.set_xlim(wmin, wmax)
#ax12.set_ylim(3e-18, 6e-16)
ax12.set_ylim(ymin/2, ymax*2)
ax12.set_yscale("log")
ax12.set_title("Extracted Spectra")

ax12.text(0.05, 0.10, "DET%s" % (det), transform=ax12.transAxes )
ax12.text(0.05, 0.05, "Gain = %.2f" % (gain), transform=ax12.transAxes )

ax13.set_xlim(wmin, wmax)
ax13.set_ylim(-0.5,4)
ax13.set_title("Input Galaxy Model/Extracted spectra [ratio]")


cb = fig1.colorbar(im, ax=axes1)
cb.set_label('Magnitude', rotation=270)
cb.ax.invert_yaxis() 

print(ymin,ymax)

#fig1.suptitle('%s DET%s Frame%i' % (det, frame), fontsize=16)
fig1.suptitle(os.path.basename(slitless_file))


plt.show()
