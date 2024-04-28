 # export PYSYN_CDBS=/Users/gwalth/data/JWST/grp/redcat/trds/

import os
print(os.environ['PYSYN_CDBS'])
import numpy as np
import pysynphot as S

from scipy import interpolate, integrate

import matplotlib.pyplot as plt

# constants
c_cm = 2.9979E10       # cm/s
h = 6.626068E-27    # cm^2*g/s
kb = 1.3806503E-16   # cm^2*g/(s^2*K) 

Ang = 1E-8          # cm
mu = 1E-4           # cm


def mag_fn(y, C):
    return -2.5*np.log10(y) - C

def ab_to_fnu(mab=None):
    """Convert AB magnitude to Fnu [erg/s/cm^2/Hz]"""
    # mab = -2.5*log10(fnu) - 48.6
    fnu = 10 ** (-0.4 * (mab + 48.6))
    return fnu
    
def fnu_to_flam(fnu=None, wav=None):
    """Convert Fnu [erg/s/cm^2/Hz] to Flam [erg/s/cm^2/Ang]"""
    flam = fnu * c_cm / (wav * Ang) ** 2 * Ang
    return flam

def flam_to_fnu(flam=None, wav=None):
    """Convert Flam [erg/s/cm^2/Ang to Fnu [erg/s/cm^2/Hz]]"""
    fnu = flam * (wav * Ang) ** 2 / c_cm / Ang
    return fnu


def calc_pivot(tran, wav):

    numerator = integrate.trapz(tran * wav, wav)
    denominator = integrate.trapz(tran / wav, wav)
    pivot = np.sqrt(numerator/denominator)

    return pivot

def calc_spectra_mag(fnu, wav, filt_trans, filt_wav, unit="Jy", area=0.785, plot=1):
    """
    fnu: np.array
    wav: np.array
    area (default) = pi * r**2 (r=0.5 arcsec)
    """
    ############
    # sed info #
    ############

    Consts = [48.6, -23.9, -31.4]
    units  = ["Jy", "uJy", "nJy"]


    # interpolate source spectrum
    S = interpolate.interp1d(wav, fnu)


    #print(filt_wav.shape)
    #print(filt_trans.shape)

    if plot:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(filt_wav, filt_trans, c="tab:orange")

        #ax2 = ax1.twinx()
        #tmp_flam = fnu_to_flam(fnu=fnu, wav=model_wav)

        #ax2.plot(wav, tmp_flam, c="tab:blue")

        plt.show()

    # interpolate filter
    F = interpolate.interp1d(filt_wav, filt_trans)

    filt_all_int = integrate.trapz(filt_trans, filt_wav)
    dw = wav[1]-wav[0]
    #filt_all_sum  = np.sum(tran)*dw
    print(filt_all_int)
    #print filt_all_sum

    # true_flux = flux_measured_in_filter * C 
    #print(wav)
    #print(filt_wav)

    i0 = np.argmin(np.abs(wav[0] - filt_wav))
    i1 = np.argmin(np.abs(wav[-1] - filt_wav))

    F_all = np.array([S(w)*F(w) for w in filt_wav[i0+1:i1-1]])

    flux_int = integrate.trapz(F_all, filt_wav[i0+1:i1-1])
    #flux_sum = np.sum(F_all)*dw

    print(flux_int)
    #print flux_sum
    flux_corr = (flux_int/filt_all_int) # Jy
    #flux_corr = (flux_int/filt_all_int)*10**-9*Jy  # nJy
    #print flux_corr           # erg/s/cm^2/Hz
    #print flux_corr/Jy        # Jy
    #print flux_corr*10**9/Jy  # nJy
    #print flux_sum/filt_all_int
    #mAB = -2.5*np.log10(flux_corr) - 48.6
    #mAB_5sig = -2.5*np.log10(flux_corr/5.) - 48.6
    #mAB_10sig = -2.5*np.log10(flux_corr/10.) - 48.6
    #mAB_area = -2.5*np.log10(flux_corr/area) - 48.6

    ind = units.index(unit)
    Const = Consts[ind]

    mAB = mag_fn(flux_corr, Const)
    print("mAB =", mAB)
    mAB_5sig = mag_fn(flux_corr/5., Const)
    print("mAB (5 sig) =", mAB_5sig)
    mAB_10sig = mag_fn(flux_corr/10., Const)
    print("mAB (10 sig) =", mAB_10sig)
    mAB_area = mag_fn(flux_corr/area, Const)
    print('mAB (mag/") =', mAB_area)
    print()

    return mAB, flux_corr

def renormalize(spec_obj, bp_obj, norm_abmag):
    fnu = flam_to_fnu(flam=spec_obj.flux, wav=spec_obj.wave)

    mAB, fnu_band = calc_spectra_mag( fnu, spec_obj.wave, bp_obj.throughput, bp_obj.wave, plot=0)
    pivot_wav = calc_pivot(bp_obj.throughput, bp_obj.wave)
    
    print(mAB, fnu_band)
    print(pivot_wav)
    
    flambda_norm = fnu_band*(c_cm/(pivot_wav*Ang)**2)*Ang  # erg/s/cm^2/Ang
    print(flambda_norm)
    
    # mab = -2.5*log10(fnu) - 48.6
    fnu = 10**((norm_abmag + 48.6)/-2.5) # erg/s/cm^2/Hz
    print(fnu)
    flambda = fnu*(c_cm/(pivot_wav*Ang)**2)*Ang  # erg/s/cm^2/Ang
    print(flambda)
    
    fnu = flam_to_fnu(flam=spec_obj.flux * flambda / flambda_norm, wav=spec_obj.wave)
    
    mAB, fnu_band = calc_spectra_mag(fnu, spec_obj.wave, bp_obj.throughput, bp_obj.wave, plot=0)
    pivot_wav = calc_pivot(bp_obj.throughput, bp_obj.wave)
    print(mAB, fnu_band)
    print(pivot_wav)
    print()
    flam_band = fnu_to_flam(fnu=fnu_band, wav=pivot_wav)
 
    return flambda, flambda_norm, pivot_wav, flam_band

#############################################################################################



vega_file = os.path.join(os.environ['PYSYN_CDBS'], 'calspec', 'alpha_lyr_stis_005.fits')
vega = S.FileSpectrum(vega_file)
vega.convert('flam')

print(vega(5000))  # photlam (internal unit)
print(vega.sample(5000))



fig = plt.figure()

ax = fig.add_subplot(111)
ax.plot(vega.wave, vega.flux)
ax.set(
    xlim=(0,12000), 
    xlabel=vega.waveunits, 
    ylabel=vega.fluxunits, 
    title=os.path.basename(vega.name)
)
plt.show()




bp1 = S.ObsBandpass('acs,wfc1,f814w')
#bp1.showfiles()

bp2 = S.ObsBandpass('wfc3,ir,f125w')
#bp2.showfiles()

bp3 = S.ObsBandpass('wfc3,ir,f160w')
#bp3.showfiles()



fig = plt.figure()

ax = fig.add_subplot(111)
ax.plot(bp1.wave, bp1.throughput)
ax.plot(bp2.wave, bp2.throughput)
ax.plot(bp3.wave, bp3.throughput)
#ax.set(
#    xlim=(0,12000), 
#    xlabel=vega.waveunits, 
#    ylabel=vega.fluxunits, 
#    title=os.path.basename(vega.name)
#)
plt.show()


#sp2 = vega.renorm(1.e-16, 'flam', bp1)

norm_abmag = 15.7
sp3 = vega.renorm(norm_abmag, 'abmag', bp3)


###########################################################################
# my method
###########################################################################

flambda_f160w, flambda_norm_f160w, pivot_wav_f160w, flam_band_f160w = renormalize(vega, bp3, norm_abmag) # F160W

fnu = flam_to_fnu(flam=vega.flux * flambda_f160w / flambda_norm_f160w, wav=vega.wave)
mAB_f814w, _ = calc_spectra_mag(fnu, vega.wave, bp1.throughput, bp1.wave, plot=0)
mAB_f125w, _ = calc_spectra_mag(fnu, vega.wave, bp2.throughput, bp2.wave, plot=0)
print(mAB_f814w)
print(mAB_f125w)

flambda_f814w, flambda_norm_f814w, pivot_wav_f814w, flam_band_f814w = renormalize(vega, bp1, mAB_f814w) # F814W
flambda_f125w, flambda_norm_f125w, pivot_wav_f125w, flam_band_f125w = renormalize(vega, bp2, mAB_f125w) # F125W

########################


fig = plt.figure()

ax = fig.add_subplot(111)

ax.semilogy(vega.wave, vega.flux, 'k', label='Vega')
ax.semilogy(sp3.wave, sp3.flux, 'g', alpha=0.4, label='pysynphot F160W renorm')
ax.semilogy(vega.wave, vega.flux* flambda_f814w / flambda_norm_f814w, 'c', alpha=0.4, label='GLW F814W renorm')
ax.semilogy(vega.wave, vega.flux* flambda_f125w / flambda_norm_f125w, 'y', alpha=0.4, label='GLW F125W renorm')
ax.semilogy(vega.wave, vega.flux* flambda_f160w / flambda_norm_f160w, 'm', alpha=0.4, label='GLW F160W renorm')
ax.scatter(
    [pivot_wav_f814w, pivot_wav_f125w, pivot_wav_f160w], 
    [flam_band_f814w, flam_band_f125w, flam_band_f160w],
    c="r",
)

ax2 = ax.twinx()
ax2.plot(bp1.wave, bp1.throughput, label="F814W", alpha=0.4)
ax2.plot(bp2.wave, bp2.throughput, label="F125W", alpha=0.4)
ax2.plot(bp3.wave, bp3.throughput, label="F160W", alpha=0.4)
ax2.legend(loc='lower right')
ax2.set_ylim(0.0, 0.6)


ax.set(
    xlim=(1100,18000), 
    xlabel=vega.waveunits, 
    ylabel=vega.fluxunits, 
)
ax.legend(loc='lower left')

plt.show()




obs = S.Observation(sp3, bp3)
obs.convert('counts')





fig = plt.figure()

ax = fig.add_subplot(111)
ax.plot(obs.wave, obs.flux, marker='x', label='native')
ax.plot(obs.binwave, obs.binflux, drawstyle='steps-mid', label='binned')
ax.set(xlim=(13000, 18000), xlabel=obs.waveunits, ylabel=obs.fluxunits,)
ax.legend(loc='best')
plt.show()




st_zpt = -2.5 * np.log10(bp3.unit_response()) - 21.1
print('STmag zeropoint for {0} is {1:.5f}'.format(bp3.name, st_zpt))
