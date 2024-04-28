# export PYSYN_CDBS=/Users/gwalth/data/JWST/grp/redcat/trds/

import os
print(os.environ['PYSYN_CDBS'])
import numpy as np
import pysynphot as S

import matplotlib.pyplot as plt

# https://www.hamamatsu.com/us/en/resources/interactive-tools/photon-flux-to-radiant-flux-calculator.html

####################################
# constants
####################################
Jy = 1.0E-23        # erg/s/cm^2/Hz
Ang = 1E-8          # cm
nm = 1E-7           # cm
c = 2.9979E10       # cm/s
#h = 6.626068E-34    # J * s
h = 6.626068E-27    # erg * s 


# https://www.hamamatsu.com/us/en/resources/interactive-tools/photon-flux-to-radiant-flux-calculator.html

##################################################################################################################################

def calc_count_rate(ABmag=0, d=2.4, wav_width=850, wav_cen=5450):
    print()
    print(f" ABmag = {ABmag}\n", f"wav_cen = {wav_cen} Ang\n", f"wav_width = {wav_width} Ang\n", f"diameter = {d} m\n")
    collect_area = np.pi*(d*100./2)**2 # cm**2
    print(collect_area,"cm^2")

    #################################################################
    fnu = 10**(-0.4*(ABmag + 48.60))                 # erg/s/cm^2/Hz
    print("Calculated from magnitude")
    print(fnu,"erg/s/cm^2/Hz")
    #flambda = fnu*Ang/((lam_obs*mu)**2/c)
    flambda = fnu*Ang/((wav_cen*Ang)**2/c)
    print(flambda,"erg/s/cm^2/Ang")
    print(flambda*wav_width,"erg/s/cm^2")
    E_phot = (h*c)/(wav_cen*Ang) # erg
    print(E_phot,"erg")
    print(flambda*wav_width/E_phot,"photons/s/cm^2")
    print(flambda/E_phot,"photons/s/cm^2/Ang")
    print()
    print(flambda*wav_width*collect_area,"erg/s")
    print(flambda*wav_width/E_phot*collect_area,"photons/s")
    print(flambda/E_phot*collect_area,"photons/s/Ang")
    print()

# expectation for Vega in Vband (Vega) = 0.0, wav_width=850 Ang, wav_cen=5450 Ang ~1000 photons/s/cm**2/Ang 

# https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html
#              U     B     V     R     I      J      H      K
#wav_cen =  [3571, 4344, 5456, 6442, 7994, 12355, 16458, 21603]
wav_cen =  [3600, 4380, 5450, 6410, 7980, 12200, 16300, 21900]
wav_width = [600,  900,  850, 1500, 1500,  2600,  2900,  4100]
AB_mag =  [0.79, -0.09, 0.02, 0.21, 0.45,  0.91,  1.39,  1.85]
#          756.1 1392.6 995.5 702.0 452.0  193.1   93.3   43.6 	

print("#"*80)
print("Compare to Paul Martini's calcs")
print("#"*80)
print()
for l, dl, ab in zip(wav_cen, wav_width, AB_mag):
    calc_count_rate(ABmag=ab, d=2.4, wav_cen=l, wav_width=dl)  # Vega to AB conversion 0.02


#####################################################################################################
# comparison calculation
#####################################################################################################

print("#"*80)
print("Compare to pysynphot")
print("#"*80)
print()
vega_file = os.path.join(os.environ['PYSYN_CDBS'], 'calspec', 'alpha_lyr_stis_005.fits')
vega = S.FileSpectrum(vega_file)
vega.convert('flam')

bp1 = S.ObsBandpass('wfc3,ir,f160w')
sp3 = vega.renorm(15.7, 'abmag', bp1)

obs = S.Observation(sp3, bp1)

print(obs.primary_area,"cm**2") # cm**2
print(obs.countrate(),"counts/s")  # counts / s 
print(obs.countrate() / obs.primary_area,"counts/s/cm**2")  # counts / s / cm**2

print(obs.effstim('flam'))

print(obs.efflam())

#######################################################################
#d = 1.2 # m
d = 2.4 # m
ABmag = 15.7
#wav_cen = 17714. # Ang
wav_cen = obs.efflam()
#wav_width = 4999. # Ang
#wav_width = bp1.photbw()
wav_width = bp1.equivwidth()
calc_count_rate(ABmag=ABmag, d=d, wav_width=wav_width, wav_cen=wav_cen)
#######################################################################

#######################################################################
# checking Euclid source
d = 1.2 # m
ABmag = 15.7
wav_cen = 17714. # Ang
wav_width = 4999. # Ang
calc_count_rate(ABmag=ABmag, d=d, wav_width=wav_width, wav_cen=wav_cen)
#######################################################################

f = 5951.
zp = 25.2 # NISP H
m = -2.5*np.log10(f) + zp
print(m)
