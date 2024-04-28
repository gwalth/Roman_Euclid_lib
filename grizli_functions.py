"""
Grizli Functions
"""
import gc
import shutil
import glob,os,sys
import pickle
from math import cos, sin, atan2, pi
import time
import textwrap

import multiprocessing as mp

import astropy.io.fits as pyfits
from astropy import wcs
from astropy.table import Table, unique, join

import numpy as np

from scipy import integrate

import matplotlib.pyplot as plt

#import grizli
from grizli import multifit
#from grizli import utils, fitting
#from grizli.pipeline import auto_script

#from pympler import tracker
#from pympler import summary
#from pympler import muppy


####################################
# constants
####################################
Jy = 1.0E-23        # erg/s/cm^2/Hz
mJy = 1e-3          # Jy
uJy = 1e-6          # Jy
Mpc = 3.086e24      # cm
Ang = 1E-8          # cm
mu = 1E-4           # cm

c_km = 2.9979E5     # km/s
c = 2.9979E10       # cm/s
h = 6.626068E-27    # cm^2*g/s           # erg s 
k = 1.3806503E-16   # cm^2*g/(s^2*K)


# https://mpdaf.readthedocs.io/en/latest/_modules/mpdaf/obj/coords.html#WCS
def rotate(wcs_obj, theta):
        """Rotate WCS coordinates to new orientation given by theta.

        Analog to ``astropy.wcs.WCS.rotateCD``, which is deprecated since
        version 1.3 (see https://github.com/astropy/astropy/issues/5175).

        Parameters
        ----------
        theta : float
            Rotation in degree.

        """
        theta = np.deg2rad(theta)
        sinq = np.sin(theta)
        cosq = np.cos(theta)
        mrot = np.array([[cosq, -sinq],
                         [sinq, cosq]])

        if wcs_obj.wcs.has_cd():    # CD matrix
            newcd = np.dot(mrot, wcs_obj.wcs.cd)
            wcs_obj.wcs.cd = newcd
            wcs_obj.wcs.set()
        elif wcs_obj.wcs.has_pc():      # PC matrix + CDELT
            newpc = np.dot(mrot, wcs_obj.wcs.get_pc())
            wcs_obj.wcs.pc = newpc
            wcs_obj.wcs.set()
        else:
            raise TypeError("Unsupported wcs type (need CD or PC matrix)")

        return wcs_obj


def calc_photflam(ZP, photplam):
        """
        ZP = (-2.5*np.log10(model.photflam_list[fi]) - 21.10 -
                          5*np.log10(model.photplam_list[fi]) + 18.6921)

        ZP = -2.5*np.log10(photflam) - 21.10 - 5*np.log10(photplam) + 18.6921
        """

        photflam = 10**(-0.4*(ZP + 21.10 - 18.6921 + 5*np.log10(photplam)))

        return photflam

def calc_ZP(photflam, photplam):

    ZP = -2.5*np.log10(photflam) - 21.10 - 5*np.log10(photplam) + 18.6921

    return ZP


def check_sims2(sim, mag_limit):

    mask = sim.catalog['MAG_AUTO'] < mag_limit

    all_cat = sim.catalog
    magcut_cat = sim.catalog[mask]

    new_ids = []
    for id in sim.object_dispersers:
        is_cgs, spectrum_1d, beam = sim.object_dispersers[id]
        if len(beam) > 0:
            new_ids.append(id)

    new_ids = np.array(new_ids)
    #print(new_ids)

    # catalog accepts either indices or boolean mask
    indices = [list(sim.catalog["NUMBER"]).index(i) for i in new_ids]
    extract_cat = sim.catalog[indices]

    print("All sources   =",len(all_cat))
    print("Magnitude cut =",len(magcut_cat))
    print("Simulated     =",len(new_ids))
    print()
    #print(len(extract_cat))

    return all_cat,magcut_cat,extract_cat

def fake_euclid_ref(tbl, ra_cen = 0.0, dec_cen = 0.0, pixel_scale = 0.3, naxis1 = 10000, 
                    naxis2 = 10000, flux_key="TU_FNU_H_NISP_MAG", mag_key="NIR_H", id_key='SOURCE_ID',
                    gain=1.0, verb=0,
                    collect_area = 9926., wav_cen = 17714., wav_width = 4999., background=1.12, 
                    exptime=100., nexp=1, readnoise=6., eff_tot = 0.8, output = "Euclid_direct.fits",
                    filt = "NISP_H", instr = "NISP", theta=0.0, photflam=1e-20):

    
    
    data = np.zeros((naxis2,naxis1)) 
    err = np.zeros((naxis2,naxis1))
    dq = np.zeros((naxis2,naxis1))
    
    w = build_wcs(naxis1 = naxis1, naxis2 = naxis2, ra_cen = ra_cen, dec_cen = dec_cen, 
                  theta=theta, pixel_scale = pixel_scale)
    
    #wcs_header = w.to_header()
    
    fluxes = []
    mags = []
    coord_offsets = []

    skipped = 0
    simulated = 0
    modified = 0
    #print(w.wcs)    
    #print(w)
    for k,row in enumerate(tbl):
        
        thumb_norm = row["THUMB"]
        src_id = row[id_key]
                
        yl,xl = thumb_norm.shape

        #world = [[tbl['RA'][k], tbl['DEC'][k]]]        
        world = [[tbl['RA_MAG'][k], tbl['DEC_MAG'][k]]]        
        pixcrd = w.wcs_world2pix(world, 0)    
        #print(world)
        #print(pixcrd)
        
        xpix, ypix = pixcrd[0]





        if xpix > 0:
            xsign = 1
        else:
            xsign = -1
    
        if ypix > 0:
            ysign = 1
        else:
            ysign = -1
    
    
        #print(xpix, ypix, xpix_new, ypix_new)
    
        # either 0.0 or 0.5 fractions of pixel
        if xl % 2:
            # if remainder
            xpix_new = int(xpix) + xsign*0.5
    
        else:
            # if no remainder
            xpix_new = int(xpix + xsign*0.5)
    
        if yl % 2:
            # if remainder
            ypix_new = int(ypix) + ysign*0.5
    
        else:
            # if no remainder
            ypix_new = int(ypix + ysign*0.5)
    
        box_hw_x = xl/2
        box_hw_y = yl/2
    
        y0 = int(ypix_new - box_hw_y)
        y1 = int(ypix_new + box_hw_y)
        x0 = int(xpix_new - box_hw_x)
        x1 = int(xpix_new + box_hw_x)
    
        dx = xpix - xpix_new
        dy = ypix - ypix_new
    
        #print(x0, x1, y0, y1, xl, yl, x1-x0, y1-y0, dx, dy)
        #print(xpix, ypix, xpix_new, ypix_new)
        


        if xpix < 0 and x0 < 0 and x1 < 0: 
            print()
            print("SKIPING object %i, outside of simulated range" % k)
            print(xpix,ypix)
            print(x0,x1,y0,y1)
            print()
            skipped += 1
            continue

        if xpix > naxis1 and x0 > naxis1 and x1 > naxis1:
            print()
            print("SKIPING object %i, outside of simulated range" % k)
            print(xpix,ypix)
            print(x0,x1,y0,y1)
            print()
            skipped += 1
            continue
              
        if ypix < 0 and y0 < 0 and y1 < 0: 
            print()
            print("SKIPING object %i, outside of simulated range" % k)
            print(xpix,ypix)
            print(x0,x1,y0,y1)
            print()
            skipped += 1
            continue

        if ypix > naxis2 and y0 > naxis2 and y1 > naxis2:
            print()
            print("SKIPING object %i, outside of simulated range" % k)
            print(xpix,ypix)
            print(x0,x1,y0,y1)
            print()
            skipped += 1
            continue

        x0t = 0
        x1t = xl
        y0t = 0
        y1t = yl

        if x0 < 0:
            print()
            print("MODIFYING object %i, edge of simulated range" % k)
            print(xpix,ypix)
            print(x0,x1,y0,y1)
            print(x0t,x1t,y0t,y1t)
            x0t = np.abs(x0)
            x0 = 0
            print(x0,x1,y0,y1)
            print(x0t,x1t,y0t,y1t)
            print()
            modified += 1
        if x1 > naxis1:
            print()
            print("MODIFYING object %i, edge of simulated range" % k)
            print(xpix,ypix)
            print(x0,x1,y0,y1)
            print(x0t,x1t,y0t,y1t)
            x1t = naxis1 - x1
            x1 = naxis1
            print(x0,x1,y0,y1)
            print(x0t,x1t,y0t,y1t)
            print()
            modified += 1
        if y0 < 0:
            print()
            print("MODIFYING object %i, edge of simulated range" % k)
            print(xpix,ypix)
            print(x0,x1,y0,y1)
            print(x0t,x1t,y0t,y1t)
            y0t = np.abs(y0)
            y0 = 0
            print(x0,x1,y0,y1)
            print(x0t,x1t,y0t,y1t)
            print()
            modified += 1
        if y1 > naxis2:
            print()
            print("MODIFYING object %i, edge of simulated range" % k)
            print(xpix,ypix)
            print(x0,x1,y0,y1)
            print(x0t,x1t,y0t,y1t)
            y1t = naxis2 - y1
            y1 = naxis2
            print(x0,x1,y0,y1)
            print(x0t,x1t,y0t,y1t)
            print()
            modified += 1
    
        fnu_Jy = tbl[flux_key][k] # Jy                
        fnu = fnu_Jy*Jy # erg/s/cm^2/Hz
        flambda = fnu*Ang/((wav_cen*Ang)**2/c) # erg/s/cm^2/Ang
        E_phot = (h*c)/(wav_cen*Ang) # erg
        flux = flambda*wav_width*collect_area*eff_tot/E_phot # photons/s
                
        fluxes.append(flux)
        mags.append(tbl[mag_key][k])
        #coord_offsets.append([src_id, dx, dy])
        
        offset_dict = {id_key: src_id, 'DX': dx, 'DY': dy}
        coord_offsets.append(offset_dict)
        flux_thumb = thumb_norm*flux
            
        if verb:
            print(k)
            #print(head['EXTNAME'])
            print(thumb_norm.shape)
            print(flux)
            print(np.sum(flux_thumb))
            print(x0,x1)
            print(y0,y1)
            print()

        
        #data[y0:y1,x0:x1] = data[y0:y1,x0:x1] + flux_thumb  # e-
        try:
            data[y0:y1,x0:x1] = data[y0:y1,x0:x1] + flux_thumb[y0t:y1t,x0t:x1t]  # e-
        #except:
        except Exception as error:
            # handle the exception
            print("An exception occurred:", error) 
            print()
            print(flux_thumb.shape)
            print(naxis1,naxis2)
            print("(xpix,ypix) =", xpix,ypix)
            print("(box_hw_x,box_hw_y) =", box_hw_x,box_hw_y)
            print("(x0,x1,y0,y1) =", x0,x1,y0,y1)
            print("(x0t,x1t,y0t,y1t) =", x0t,x1t,y0t,y1t)

            print("data cutout =", data[y0:y1,x0:x1].shape)
            print("flux_thumb cutout =", flux_thumb[y0t:y1t,x0t:x1t].shape)
            sys.exit(0)

        sys.stdout.write("Progress: %d   \r" % (k+1))
        sys.stdout.flush()

        simulated += 1

    print()

    print("Sources simulated =", simulated)
    print("Sources modified =", modified)
    print("Sources skipped =", skipped)

    #hdu = pyfits.PrimaryHDU(header=wcs_header)  
    #hdu.name = 'PRIMARY'
    hdu = pyfits.PrimaryHDU()  
    hdu.name = 'PRIMARY'
    
    cd_mat = cd_matrix(-1*pixel_scale/3600, pixel_scale/3600, theta) 
    
    wcs_dict = {'CRPIX1': naxis1/2., 'CRPIX2': naxis2/2., 'CRVAL1': ra_cen, 'CRVAL2': dec_cen,
                'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CUNIT1': 'deg', 'CUNIT2': 'deg',
                'NAXIS1': naxis1, 'NAXIS2': naxis2, 'CD1_1': cd_mat[0,0], 'CD1_2': cd_mat[0,1], 
                'CD2_1': cd_mat[1,0], 'CD2_2': cd_mat[1,1],} # 'PC1_1': 1.0, 'PC2_2': 1.0}
        
    hdu_sci = pyfits.ImageHDU(data)  
    #hdu_sci = pyfits.ImageHDU(data, header=wcs_header)  
    #hdu_sci.name = 'SCI'
    hdu_sci.name = 'REF'

    hdu_err = pyfits.ImageHDU(err)
    #hdu_err = pyfits.ImageHDU(err, header=wcs_header)
    hdu_err.name = 'ERR'
    hdu_dq = pyfits.ImageHDU(dq)
    #hdu_dq = pyfits.ImageHDU(dq, header=wcs_header)
    hdu_dq.name = 'DQ'

    [hdu_sci.header.set(wcs_key, value=wcs_dict[wcs_key]) for wcs_key in wcs_dict]
    [hdu_err.header.set(wcs_key, value=wcs_dict[wcs_key]) for wcs_key in wcs_dict]
    [hdu_dq.header.set(wcs_key, value=wcs_dict[wcs_key]) for wcs_key in wcs_dict]
 
    hdu1 = pyfits.HDUList([hdu,hdu_sci,hdu_err,hdu_dq])
    hdu1.info()

    # hdu e-  --> e-/s
    #hdu2 = add_noise(hdu1, background=background, exptime=exptime, 
    #                nexp=nexp, readnoise=readnoise)
    hdu2 = add_noise(hdu1, background=background, exptime=exptime, 
                     nexp=nexp, readnoise=readnoise, img_ext='REF')
       
    ext = 0
    hdu2[ext].header['FILTER'] = filt
    hdu2[ext].header['INSTRUME'] = instr
    hdu2[ext].header['EXPTIME'] = exptime 
    hdu2[ext].header['PHOTPLAM'] = wav_cen
    hdu2[ext].header['PHOTFLAM'] = photflam
    ##hdu[ext].header['EXTVER'] = ext

    ext = 1
    hdu2[ext].header['FILTER'] = filt
    hdu2[ext].header['INSTRUME'] = instr
    hdu2[ext].header['EXTVER'] = ext
    hdu2[ext].header['EXPTIME'] = exptime   
    hdu2[ext].header['PHOTPLAM'] = wav_cen
    hdu2[ext].header['PHOTFLAM'] = photflam
    
    hdu2.writeto("../" + output, overwrite=True)
            
            
    if verb:
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111)
        ax.hist(np.log10(fluxes),bins=30)
        ax.set_yscale("log")
        plt.show()
            
    return fluxes, mags, coord_offsets

def print_header(w):
    wcs_header = w.to_header()
    for l in textwrap.wrap(str(wcs_header), width=80): print(l)

def build_wcs(naxis1 = 10000, naxis2 = 10000, ra_cen = 0.0, dec_cen = 0.0, 
              pixel_scale = 0.3, theta = 0.0):
    # Create a new WCS object.  The number of axes must be set
    # from the start
    
    cd_mat = cd_matrix(-1*pixel_scale/3600, pixel_scale/3600, theta) 
    
    
    wcs_dict = {'CRPIX1': naxis1/2., 'CRPIX2': naxis2/2., 'CRVAL1': ra_cen, 'CRVAL2': dec_cen,
                'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CUNIT1': 'deg', 'CUNIT2': 'deg',
                'NAXIS1': naxis1, 'NAXIS2': naxis2, 'CD1_1': cd_mat[0,0], 'CD1_2': cd_mat[0,1], 
                'CD2_1': cd_mat[1,0], 'CD2_2': cd_mat[1,1],} # 'PC1_1': 1.0, 'PC2_2': 1.0}
    
    w = wcs.WCS(wcs_dict)    
    w.wcs.set()
    
    #print_header(w)
    #w.array_shape = [naxis1, naxis2]
    
    #w.to_header()
    
    #w.printwcs()
    #print()
    
    #print(w.to_header_string())
    #print(wcs_header)
    #hdu = w.to_fits()
    
    return w

def cd_matrix(cdelt1, cdelt2, theta):
    cd1_1 =  cdelt1 * np.cos(theta*np.pi/180.)
    cd1_2 = -cdelt2 * np.sin(theta*np.pi/180.)
    cd2_1 =  cdelt1 * np.sin(theta*np.pi/180.)
    cd2_2 =  cdelt2 * np.cos(theta*np.pi/180.)

    return np.array([[cd1_1, cd1_2],
                     [cd2_1, cd2_2]])

def total_image(files, output="total.fits", method="sum", img_ext='SCI'):
    # assuming images are all the same size
    
    all_sci = []
    all_err = []
    all_dq  = []
    
    for f in files:
        
        hdu = pyfits.open(f)
        sci = hdu[img_ext].data
        err = hdu['ERR'].data
        dq  = hdu['DQ'].data
        
        head = hdu[img_ext].header
        
        all_sci.append(sci)
        all_err.append(err)
        all_dq.append(err)


    all_sci = np.array(all_sci)
    all_err = np.array(all_err)

    print(all_sci)
    print(all_sci.shape)
    
    N = len(all_sci)
        
    if method == "sum":
        total_sci = np.sum(all_sci, axis=0)
        total_err = np.sqrt(np.sum(all_err**2, axis=0))
        total_dq = all_dq[0] # kludge - need to think about proper dq
        
    elif method == "mean":
        total_sci = np.mean(all_sci, axis=0)
        total_err = np.sqrt(np.sum(all_err**2, axis=0))/N
        total_dq = all_dq[0] # kludge - need to think about proper dq

        
    hdu = pyfits.PrimaryHDU(header=head)
    hdu.name = 'PRIMARY'

    hdu_sci = pyfits.ImageHDU(total_sci, header=head)
    hdu_sci.name = img_ext
    hdu_err = pyfits.ImageHDU(total_err, header=head)
    hdu_err.name = 'ERR'
    hdu_dq = pyfits.ImageHDU(total_dq, header=head)
    hdu_dq.name = 'DQ'

    hdu1 = pyfits.HDUList([hdu,hdu_sci,hdu_err,hdu_dq])
    hdu1.info()
    hdu1.writeto(output, overwrite=True)

    #print(files)

def calc_pivot(wav,tran):
    
    numerator = integrate.trapz(tran * wav, wav)
    denominator = integrate.trapz(tran / wav, wav)    
    pivot = np.sqrt(numerator/denominator)
    
    return pivot

def euclid_det():

    all_det = []

    for i in range(0,4): 
        for j in range(0,4):
            print(i+1,j+1)
            num = '%i%i' % (i+1,j+1)
            det = 'DET%s' % (num)       
            all_det.append(det)
    print(all_det)
    return all_det

def remove_basic_wcs(head, alt=''):
 
    # alt is for alternate WCS (i.e. A, B, C,...)   

    keys_to_delete = [
        "CRVAL1",
        "CRVAL2",
        "CRPIX1",
        "CRPIX2",
        "CTYPE1",
        "CTYPE2",
        "CUNIT1",
        "CUNIT2",
        "CD1_1",
        "CD1_2",
        "CD2_1",
        "CD2_2",
    ]

    keys_to_delete = [k + alt for k in keys_to_delete]

    keys_to_delete += ["CUNIT1", "CUNIT2"]
    
    remove_list_of_header_keys(head, keys_to_delete)


def remove_list_of_header_keys(head, keys_to_delete):
    for k in keys_to_delete:
        try:
            del head[k]
        except Exception as error:
            print("An exception occurred:", error)

def euclid_wcs(head):

    wlist = wcs.find_all_wcs(head)

    for w in wlist:

        if 'Y_MOSA' in w.wcs.ctype and 'Z_MOSA' in w.wcs.ctype:
            #w.wcs.alt = 'S'
            continue
        elif 'RA---TAN' in w.wcs.ctype and 'DEC--TAN' in w.wcs.ctype:
            w.wcs.alt = ' '

        newhead = w.to_header()
    
        newhead["CD1_1"] = newhead["PC1_1"]
        newhead["CD1_2"] = newhead["PC1_2"]
        newhead["CD2_1"] = newhead["PC2_1"]
        newhead["CD2_2"] = newhead["PC2_2"]
    
        keys_to_delete = [
            "PC1_1",
            "PC1_2",
            "PC2_1",
            "PC2_2",
            "CDELT1",
            "CDELT2",
            "LONPOLE",
            "LATPOLE",
            "MJDREF",
            "RADESYS",
        ]
    
        remove_list_of_header_keys(newhead, keys_to_delete)
    
    return newhead
    

def write_individual_slitless(slitless_file, file_str = "Euclid_DET%s_slitless.fits", rot=0):
    
    # Read in Euclid slitless multi-extension fits and write single frame for each detector

    # Loops through all detectors

    hdu_old = pyfits.open(slitless_file)
    print(hdu_old.info())

    # for now dq is ones
    dq = np.ones((2048,2048)) 

    all_slitless = []

    for i in range(0,4): 
        for j in range(0,4):
            #print(i+1,j+1)
            num = '%i%i' % (i+1,j+1)
            #det = 'DET' + num
        
            #conf_file = 'NISP-GLWv2-DET%s' % (num)
            #slitless = "Euclid_DET%s_slitless.fits"  % (num)
            slitless = file_str % (num)
            #print(conf_file)
            print(slitless)
        
            all_slitless.append(slitless)

            # test to read conf_file in grismconf.py
            #det = conf_file.split("-")[-1][-2:]
            #print(det)
        
            #head = hdu_old['DET%s.SCI' % (num)].header
            data = hdu_old['DET%s.SCI' % (num)].data
            err = hdu_old['DET%s.CHI2' % (num)].data



            # fixes WCS issue in new sims        
            head_wcs = euclid_wcs(hdu_old['DET%s.SCI' % (num)].header)
            remove_basic_wcs(hdu_old['DET%s.SCI' % (num)].header, alt='S')
            print(head_wcs)

            if rot:

                # recenter
                wcs_obj = wcs.WCS(head_wcs)
                print(wcs_obj)
                print(wcs_obj.wcs.cd)

                naxis1 = hdu_old['DET%s.SCI' % (num)].header["NAXIS1"]
                naxis2 = hdu_old['DET%s.SCI' % (num)].header["NAXIS2"]
                crpix1 = naxis1/2.
                crpix2 = naxis2/2.
                print(crpix1,crpix2)
                sky = wcs_obj.pixel_to_world(crpix1,crpix2)
                print(sky)
                crval1 = sky.ra.value
                crval2 = sky.dec.value
 
                head_wcs["CRVAL1"] = crval1
                head_wcs["CRVAL2"] = crval2
                head_wcs["CRPIX1"] = crpix1 + 1 # not sure why +1???
                head_wcs["CRPIX2"] = crpix2



                wcs_obj = wcs.WCS(head_wcs)

                new_wcs = rotate(wcs_obj, 90)
                k = 1
                #if num in ['11','12','13','14','21','22','23','24']:
                #    new_wcs = rotate(wcs_obj, 90)
                #    k = 1

                #elif num in ['31','32','33','34','41','42','43','44']:
                #    new_wcs = rotate(wcs_obj, -90)
                #    k = -1

                # rotate
                print("rotating 90 degress")
                print(num, k)
                data = np.rot90(data, k=k)
                err = np.rot90(err, k=k)


                print(new_wcs)
                print(new_wcs.wcs.cd)
                cd_mat = new_wcs.wcs.cd

                wcs_dict = {'CD1_1': cd_mat[0,0], 'CD1_2': cd_mat[0,1], 
                            'CD2_1': cd_mat[1,0], 'CD2_2': cd_mat[1,1],}

                [head_wcs.set(wcs_key, value=wcs_dict[wcs_key]) for wcs_key in wcs_dict]



            hdu_old[0].header.update(head_wcs)
            hdu_old['DET%s.SCI' % (num)].header.update(head_wcs)

            head0 = hdu_old[0].header
            head = hdu_old['DET%s.SCI' % (num)].header

            hdu = pyfits.PrimaryHDU(header=head0)  
            hdu.name = 'PRIMARY'
    
            hdu_sci = pyfits.ImageHDU(data, header=head)  
            hdu_sci.name = 'SCI'
            hdu_err = pyfits.ImageHDU(err, header=head)
            hdu_err.name = 'ERR'
            hdu_dq = pyfits.ImageHDU(dq.copy(), header=head)
            hdu_dq.name = 'DQ'

            hdu1 = pyfits.HDUList([hdu, hdu_sci, hdu_err, hdu_dq])
            hdu1.writeto(slitless, overwrite=True)
            hdu1.info()

            hdu1.close()

            #hdu[ext].header['INSTRUME'] = 'NISP-DET%i%i-GLWv2' % (i+1,j+1)

    hdu_old.close()
            
    return all_slitless

def plot_slitless(slitless_file, verb=0, vmin=200, vmax=3000):
    hdu = pyfits.open(slitless_file)
    if verb:
        print(hdu[0].header)
        hdu.info()

    fig = plt.figure(figsize=(8,8))

    k = 0
    for i in range(0,4): 
        for j in range(0,4):
            #print(i+1,j+1)
            num = '%i%i' % (i+1,j+1)
            name = 'DET%s.SCI' % (num)
            #name = 'DET%i%i.CHI2' % (i+1,j+1)
        
            img = hdu[name].data
            head = hdu[name].header
            if verb:
                print(img.shape)
                print(name)
                print(head['CRVAL1'])
                print(head['CRVAL2'])
                print(head['CRPIX1'])
                print(head['CRPIX2'])
        
                print(np.min(img))
                print(np.max(img))
                print(np.mean(img))
                print(np.median(img))
                print(np.std(img))
                print()
                #print(hdu[name].header)
                #print()
        
            p = fig.add_subplot(4,4,k+1)
            p.imshow(img,vmin=vmin,vmax=vmax)
        
            p.text(0.05,0.9,name,transform=p.transAxes)
            #p.hist(img.flatten(), bins=100)
            #p.set_yscale("log")
        
            k += 1
        
    plt.show()

def read_slitless_headers(slitless_file, verb=0, plot=0):

    hdu = pyfits.open(slitless_file)
    if verb:
        print(hdu[0].header)
        hdu.info()

    spec_exptime = hdu[0].header['exptime']

    if plot:
        fig = plt.figure(figsize=(8,8))
        ax1 = fig.add_subplot(111)

    heads = {}

    for i in range(4):
        for j in range(4):
            name = 'DET%i%i.SCI' % (i+1,j+1)
            #name = 'DET%i%i.CHI2' % (i+1,j+1)
        
            #img = hdu[name].data
            head = hdu[name].header
            heads[name] = head
            
            x = head['CRPIX1']
            y = head['CRPIX2']
            
            if verb:
                print(head)
                print(name)
                print(head['CRVAL1'])
                print(head['CRVAL2'])
                print(head['CRPIX1'])
                print(head['CRPIX2'])
                print(x,y)
                print()
            
            if plot:
                ax1.scatter(x,y)
                ax1.text(x,y,name)
    
    if plot:
        plt.show()
        
    return heads

def map_src_to_det(primer, heads, thumb_temp = 'NIS_catalog_file_%s.thm.beamA.fits',
                   wcs_temp = "../Euclid_DET%s_wcs.fits", plot=1):

    # Loops through all detectors

    det_dict = {}

    for i in range(0,4): 
        for j in range(0,4):
            print(i+1,j+1)
            num = '%i%i' % (i+1,j+1)

            head = heads['DET%s.SCI' % (num)]

            nhdu = pyfits.PrimaryHDU()                                                              
            nhdu.header = head  
            nhdu.header.tofile(wcs_temp % (num),overwrite=True)
            #print(head)
            #print()
            #print(nhdu)                                                                            
                                                                                                                                                                               
            w = wcs.WCS(nhdu)  

            hdu = pyfits.open(thumb_temp % (num))
            sources = [h.header['EXTNAME'] for h in hdu[1:]]
            det_thumb_tbl = Table([sources], names=['SOURCE_ID'])
        
            det_tbl = join(det_thumb_tbl, primer, keys='SOURCE_ID')

            world = [[ra,dec] for ra,dec in det_tbl['RA','DEC']]
            pixcrd = w.wcs_world2pix(world, 0)                                              

            det_tbl['X_PIX'] = pixcrd[:,0]
            det_tbl['Y_PIX'] = pixcrd[:,1]
        
            if plot:
                fig = plt.figure()
        
                p1 = fig.add_subplot(121)
                p1.scatter(det_tbl['X_PIX'],det_tbl['Y_PIX'])
                p1.set_aspect("equal")
                #p1.set_xlim()
        
                p2 = fig.add_subplot(122)
                p2.scatter(det_tbl['X_PIX'],det_tbl['Y_PIX'])
        
                x0,x1 = p2.get_xlim()
        
                p2.plot([x0,x1],[0,0],"--",c="r")
                p2.plot([x0,x1],[2048,2048],"--",c="r")        
                p2.set_aspect("equal")        
        
                plt.show()
        
            print("N =",sources)
            print()
            print()
        
            det_dict[num] = det_tbl

    print(det_dict)
    #return det_tbl, det_dict
    return det_dict

def wcs_pixel_scale(file, ext=1):

    hdu = pyfits.open(file)
    hdu.info()
    head = hdu[ext].header

    naxis2 = head["NAXIS2"]
    crpix1 = head["CRPIX1"]
    crpix2 = head["CRPIX2"]
    crval1 = head["CRVAL1"]
    crval2 = head["CRVAL2"]
    ctype1 = head["CTYPE1"]
    ctype2 = head["CTYPE2"]
    
    cd1_1 = head["CD1_1"]
    cd1_2 = head["CD1_2"]
    cd2_1 = head["CD2_1"]
    cd2_2 = head["CD2_2"]
    
    crota1 = atan2( cd2_1, cd1_1)  # radians
    crota2 = atan2(-cd1_2, cd2_2)  # radians
    cdelt1 = cd1_1/cos(crota1)     # deg/pix
    cdelt2 = cd2_2/cos(crota1)     # deg/pix

    print('cdelt1 = %.4f "/pixel' % (cdelt1*3600))
    print('cdelt2 = %.4f "/pixel' % (cdelt2*3600))
    print('crota1 = %.4f deg' % (crota1*180/pi))
    print('crota2 = %.4f deg' % (crota2*180/pi))
    print()


def create_circular_mask(h, w, center, radius):

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

def check_sims(sim, mag_limit):
    
    mask = sim.catalog['MAG_AUTO'] < mag_limit
    
    all_cat = sim.catalog
    magcut_cat = sim.catalog[mask]

    new_ids = []
    for id in sim.object_dispersers:
        is_cgs, spectrum_1d, beam = sim.object_dispersers[id]
        if len(beam) > 0: 
            new_ids.append(id)
        
    new_ids = np.array(new_ids)
    #print(new_ids)
        
    extract_cat = sim.catalog[new_ids - 1] 
    
    print("All sources   =",len(all_cat))
    print("Magnitude cut =",len(magcut_cat))
    print("Simulated     =",len(new_ids))
    print()
    #print(len(extract_cat))
    
    return all_cat,magcut_cat,extract_cat


def display_grizli_roman(root, id, w0=0.9, w1=2.5, labels=1, y0=-1, y1=-1):

    from matplotlib.gridspec import GridSpec
    from matplotlib.ticker import MultipleLocator
    
    # constants
    fontsize = 8
    #y0 = -1e-18
    #y1 = 3.5e-18
    lw = 2
    

    f_full = '{0}_{1:05d}.full.fits'.format(root, id)
    #print(f_full)
    full_hdu = pyfits.open(f_full)
    #print(full_hdu.info())
    head = full_hdu[0].header
    zfit_stack = Table(full_hdu['ZFIT_STACK'].data)
    zfit_head = full_hdu['ZFIT_STACK'].header
    templ = Table(full_hdu['TEMPL'].data)
    print()
    
    #print(head)
    ndfilts = head["NDFILT"] # number of direct image filters
    for i in range(ndfilts):
        print(head["DFILT%02d" % (i+1)])
        
        direct = full_hdu[5+(2*i)].data
        #print(direct)
        #print(direct.shape)    
    
    #head["REDSHIFT"] # redshift
    #head["NUMLINES"] # number of emission lines
    

    f_1d = '{0}_{1:05d}.1D.fits'.format(root, id)
    #print(f_1d)
    oned_hdu = pyfits.open(f_1d)
    #print(oned_hdu[1].header)
    print(oned_hdu.info())
    grism = Table(oned_hdu['GRISM'].data)
    grism.show_in_notebook()
    print(grism.colnames)
    #print()
    
    
    
    f_2d = '{0}_{1:05d}.stack.fits'.format(root, id)
    #print(f_2d)
    twod_hdu = pyfits.open(f_2d)
    #print(twod_hdu.info())
    #print()

    
    fig = plt.figure(figsize=(12,10))
    #           row column
    gs = GridSpec(3, 4, hspace=0.4)

    p1 = fig.add_subplot(gs[:1, 0:2])
    p2 = fig.add_subplot(gs[:1, 2:])
    p3 = fig.add_subplot(gs[1, :])
    p4 = fig.add_subplot(gs[2, :])


    # Z_MAP, CHIMIN and DOF
    
    p1.text(0.95, 0.96, root + '\n'+'ID={0:<5d}  z={1:.4f}'.format(id, zfit_head['z_map']), ha='right', va='top', transform=p1.transAxes, fontsize=9)

    zmi, zma = zfit_stack['zgrid'].min(), zfit_stack['zgrid'].max()
    if (zma-zmi) > 5:
        ticks = np.arange(np.ceil(zmi), np.floor(zma)+0.5, 1)
        lz = np.log(1+zfit_stack['zgrid'])
        p1.plot(lz, np.log10(zfit_stack['pdf']), color='k')
        p1.set_xticks(np.log(1+ticks))
        p1.set_xticklabels(np.cast[int](ticks))
        p1.set_xlim(lz.min(), lz.max())
    else:
        p1.plot(zfit_stack['zgrid'], np.log10(zfit_stack['pdf']), color='k')
        p1.set_xlim(zmi, zma)

    p1.set_xlabel(r'$z$')
    p1.set_ylabel(r'$\log\ p(z)$'+' / '+ r'$\chi^2=\frac{{{0:.0f}}}{{{1:d}}}={2:.2f}$'.format(zfit_head['chimin'], zfit_head['DoF'], zfit_head['chimin']/zfit_head['DoF']))
    p1.set_yticks([1,4,9,16,25])

    pzmax = np.log10(zfit_stack['pdf'].max())
    p1.set_ylim(pzmax-6, pzmax+0.9)
    p1.grid()
    p1.yaxis.set_major_locator(MultipleLocator(base=1))
    
    
    
    #pz_max = np.log10(zfit_stack['pdf'].max())
    
    #p1.plot(zfit_stack['zgrid'], np.log10(zfit_stack['pdf']), label='Stacked')
    #p1.set_xlim(0.0, 3); 
    ##p1.semilogy(); 
    #p1.grid()
    ##p1.set_ylim(1.e-50, 1e4)
    #p1.set_ylim(pz_max-6, pz_max+0.9)
    #p1.set_xlabel('z'); p1.set_ylabel('PDF(z)'); #plt.legend()

    p2.plot(zfit_stack['zgrid'], zfit_stack['risk'], label='Stacked')
    p2.set_xlim(0.0, 3); p2.semilogy(); p2.grid()
    p2.set_xlabel('z'); p2.set_ylabel('risk'); #p3.legend()
    
    #print('Continuum template, cont1d: ', templ['continuum'].__class__)
    #print(templ.colnames)
    
    
    if y0 < 0: y0 = np.min(templ['continuum'] + templ['full'])
    if y1 < 0: y1 = np.max(templ['continuum'] + templ['full'])
    
    
    p3.plot(templ['wave']/1.e4, templ['continuum'], label='continuum')
    p3.plot(templ['wave']/1.e4, templ['full'], label='total')
    p3.set_xlim(w0, w1); p3.set_ylim(y0,y1);#p3.set_ylim(0,1.e-17) 
    p3.grid()
    p3.set_xlabel(r'$\lambda$ (microns)')
    p3.set_ylabel(r'F$_\lambda$ (erg/s/cm$^2$/$\AA$)')
    #p3.legend()
    
    
        
    z0 = zfit_head['z_map']

    for line,wemit in emlines:
        wobs = (1+z0)*wemit
        #FWHM = vel * wobs / c_km
        if wobs/1e4 > w0 and wobs/1e4 < w1:

            p4.plot([wobs/1e4,wobs/1e4],[y0,y1],":",lw=lw,c="b")
        
            if labels:
                p4.text(wobs/1e4,0.7*y1,line,fontsize=fontsize,
                        rotation='vertical',
                        horizontalalignment='center',
                        verticalalignment='center')
    
        
        
    p4.errorbar(grism['wave']/1e4, grism['flux']/grism['flat'], yerr=grism["err"]/grism['flat'], 
                color="g", marker='.', linestyle='None', alpha=0.5) 
                
    p4.plot(grism['wave']/1e4, grism['line']/grism['flat'], color="r", alpha=0.5)  
    
    p4.plot(grism['wave']/1e4, grism['contam']/grism['flat'], color="b", alpha=0.5)
    
    p4.set_xlim(w0, w1); p4.set_ylim(y0,y1);
    #p4.semilogy(); 
    p4.grid()
    p4.set_xlabel(r'$\lambda$ (microns)')
    p4.set_ylabel(r'F$_\lambda$ (erg/s/cm$^2$/$\AA$)')
    
    # Gabe's routine
    multifit.show_drizzle_HDU(twod_hdu)

def display_grizli_jwst(root, id, w0=0.8, w1=2.1, labels=1, y0=1, y1=-1, z=None):

    from matplotlib.gridspec import GridSpec
    from matplotlib.ticker import MultipleLocator
    
    # constants
    fontsize = 8
    #y0 = -1e-18
    #y1 = 3.5e-18
    lw = 2
    

    f_full = '{0}_{1:05d}.full.fits'.format(root, id)
    #print(f_full)
    full_hdu = pyfits.open(f_full)
    #print(full_hdu.info())
    head = full_hdu[0].header
    zfit_stack = Table(full_hdu['ZFIT_STACK'].data)
    zfit_head = full_hdu['ZFIT_STACK'].header
    templ = Table(full_hdu['TEMPL'].data)
    print()
    
    #print(head)
    ndfilts = head["NDFILT"] # number of direct image filters
    for i in range(ndfilts):
        print(head["DFILT%02d" % (i+1)])
        
        direct = full_hdu[5+(2*i)].data
        #print(direct)
        #print(direct.shape)    
    
    #head["REDSHIFT"] # redshift
    #head["NUMLINES"] # number of emission lines
    

    f_1d = '{0}_{1:05d}.1D.fits'.format(root, id)
    #print(f_1d)
    oned_hdu = pyfits.open(f_1d)
    print(oned_hdu[1].header)
    print(oned_hdu.info())
    f115w = Table(oned_hdu['F115W'].data)
    f200w = Table(oned_hdu['F200W'].data)
    f115w.show_in_notebook()
    print(f115w.colnames)
    #print(f200w.colnames)
    #print()
    
    f_2d = '{0}_{1:05d}.stack.fits'.format(root, id)
    #print(f_2d)
    twod_hdu = pyfits.open(f_2d)
    #print(twod_hdu.info())
    #print()

    
    fig = plt.figure(figsize=(12,10))
    #           row column
    gs = GridSpec(3, 4, hspace=0.4)

    p1 = fig.add_subplot(gs[:1, 0:2])
    p2 = fig.add_subplot(gs[:1, 2:])
    p3 = fig.add_subplot(gs[1, :])
    p4 = fig.add_subplot(gs[2, :])


    # Z_MAP, CHIMIN and DOF
    
    p1.text(0.95, 0.96, root + '\n'+'ID={0:<5d}  z={1:.4f}'.format(id, zfit_head['z_map']), ha='right', va='top', transform=p1.transAxes, fontsize=9)

    zmi, zma = zfit_stack['zgrid'].min(), zfit_stack['zgrid'].max()
    if (zma-zmi) > 5:
        ticks = np.arange(np.ceil(zmi), np.floor(zma)+0.5, 1)
        lz = np.log(1+zfit_stack['zgrid'])
        p1.plot(lz, np.log10(zfit_stack['pdf']), color='k')
        p1.set_xticks(np.log(1+ticks))
        p1.set_xticklabels(np.cast[int](ticks))
        p1.set_xlim(lz.min(), lz.max())
    else:
        p1.plot(zfit_stack['zgrid'], np.log10(zfit_stack['pdf']), color='k')
        p1.set_xlim(zmi, zma)

    p1.set_xlabel(r'$z$')
    p1.set_ylabel(r'$\log\ p(z)$'+' / '+ r'$\chi^2=\frac{{{0:.0f}}}{{{1:d}}}={2:.2f}$'.format(zfit_head['chimin'], zfit_head['DoF'], zfit_head['chimin']/zfit_head['DoF']))
    p1.set_yticks([1,4,9,16,25])

    pzmax = np.log10(zfit_stack['pdf'].max())
    p1.set_ylim(pzmax-6, pzmax+0.9)
    p1.grid()
    p1.yaxis.set_major_locator(MultipleLocator(base=1))
    
    
    
    #pz_max = np.log10(zfit_stack['pdf'].max())
    
    #p1.plot(zfit_stack['zgrid'], np.log10(zfit_stack['pdf']), label='Stacked')
    #p1.set_xlim(0.0, 3); 
    ##p1.semilogy(); 
    #p1.grid()
    ##p1.set_ylim(1.e-50, 1e4)
    #p1.set_ylim(pz_max-6, pz_max+0.9)
    #p1.set_xlabel('z'); p1.set_ylabel('PDF(z)'); #plt.legend()

    p2.plot(zfit_stack['zgrid'], zfit_stack['risk'], label='Stacked')
    p2.set_xlim(0.0, 3); p2.semilogy(); p2.grid()
    p2.set_xlabel('z'); p2.set_ylabel('risk'); #p3.legend()
    
    #print('Continuum template, cont1d: ', templ['continuum'].__class__)
    #print(templ.colnames)
    
    
    if y0 == -1: y0 = np.min(templ['continuum'] + templ['full'])
    if y1 == -1: y1 = np.max(templ['continuum'] + templ['full'])
    
    
    p3.plot(templ['wave']/1.e4, templ['continuum'], label='continuum')
    p3.plot(templ['wave']/1.e4, templ['full'], label='total')
    p3.set_xlim(w0, w1); p3.set_ylim(y0,y1);#p3.set_ylim(0,1.e-17) 
    p3.grid()
    p3.set_xlabel(r'$\lambda$ (microns)')
    p3.set_ylabel(r'F$_\lambda$ (erg/s/cm$^2$/$\AA$)')
    #p3.legend()
    
    
    if z: z0 = z
    else: z0 = zfit_head['z_map']
    p4.text(0.95, 0.96, 'z={0:.4f}'.format(z), ha='right', va='top', transform=p4.transAxes, fontsize=9)

    for line,wemit in emlines:
        wobs = (1+z0)*wemit
        #FWHM = vel * wobs / c_km
        if wobs/1e4 > w0 and wobs/1e4 < w1:

            p4.plot([wobs/1e4,wobs/1e4],[y0,y1],":",lw=lw,c="b")
        
            if labels:
                p4.text(wobs/1e4,0.7*y1,line,fontsize=fontsize,
                        rotation='vertical',
                        horizontalalignment='center',
                        verticalalignment='center')
    
        
    p4.errorbar(f115w['wave']/1e4, f115w['flux']/f115w['flat'], yerr=f115w["err"]/f115w['flat'], 
                color="g", marker='.', linestyle='None', alpha=0.5)
    p4.errorbar(f200w['wave']/1e4, f200w['flux']/f200w['flat'], yerr=f200w["err"]/f200w['flat'], 
                color="orange", marker='.', linestyle='None', alpha=0.5)
                
    p4.plot(f115w['wave']/1e4, f115w['line']/f115w['flat'], color="r", alpha=0.5)
    p4.plot(f200w['wave']/1e4, f200w['line']/f200w['flat'], color="r", alpha=0.5)
    
    p4.plot(f115w['wave']/1e4, f115w['contam']/f115w['flat'], color="b", alpha=0.5)
    p4.plot(f200w['wave']/1e4, f200w['contam']/f200w['flat'], color="b", alpha=0.5)
    
    p4.set_xlim(w0, w1); p4.set_ylim(y0,y1);
    #p4.semilogy(); 
    p4.grid()
    p4.set_xlabel(r'$\lambda$ (microns)')
    p4.set_ylabel(r'F$_\lambda$ (erg/s/cm$^2$/$\AA$)')
    
    # Gabe's routine
    multifit.show_drizzle_HDU(twod_hdu)


# dispersers=["GR"]
# dispersers=["GRISM"]
# dispersers=["F115W","F150W","F200W"]
#def display_grizli(root, id, w0=0.8, w1=1.7, labels=1, y0=-1, y1=-1, z_in=0, path="", lw=2, 
#                   fontsize=8, dispersers=['G102','G141'], yscale=1.2):

#prefix = '{0}_{1:05d}'.format(root, id)

def display_grizli(prefix, w0=0.77, w1=1.73, labels=1, y0=-1, y1=-1, z_in=0, path="", lw=2, 
                   fontsize=8, dispersers=['G102','G141'], yscale=1.2, fig=None, cmap='viridis_r',
                   scale_size=1):

    from matplotlib.gridspec import GridSpec
    from matplotlib.ticker import MultipleLocator
    import matplotlib.colors as mcolors
    
    colors = list(mcolors.TABLEAU_COLORS)

    #f_full = '{0}_{1:05d}.full.fits'.format(root, id)
    f_full = os.path.join(path, prefix + '.full.fits')
    #print(f_full)
    full_hdu = pyfits.open(f_full)
    #full_hdu = pyfits.open(os.path.join(path,f_full))
    #print(full_hdu.info())
    head = full_hdu[0].header
    zfit_stack = Table(full_hdu['ZFIT_STACK'].data)
    zfit_head = full_hdu['ZFIT_STACK'].header
    templ = Table(full_hdu['TEMPL'].data)
    print()
    
    #print(head)
    ndfilts = head["NDFILT"] # number of direct image filters
    for i in range(ndfilts):
        print(head["DFILT%02d" % (i+1)])
        
        direct = full_hdu[5+(2*i)].data
        #print(direct)
        #print(direct.shape)    
    
    #head["REDSHIFT"] # redshift
    #head["NUMLINES"] # number of emission lines
    

    #f_1d = '{0}_{1:05d}.1D.fits'.format(root, id)
    f_1d = os.path.join(path, prefix + '.1D.fits')
    #print(f_1d)
    oned_hdu = pyfits.open(f_1d)
    #oned_hdu = pyfits.open(os.path.join(path,f_1d))
    print(oned_hdu[1].header)
    print(oned_hdu.info())


    # load 1D data
    tbl = {}
    for i,disp in enumerate(dispersers):
        tbl[i] = Table(oned_hdu[disp].data)

    #g102 = Table(oned_hdu['G102'].data)
    #g141 = Table(oned_hdu['G141'].data)
    #g102.show_in_notebook()
    #print(g102.colnames)
    #print(g141.colnames)
    #print()
    
    #f_2d = '{0}_{1:05d}.stack.fits'.format(root, id)
    f_2d = os.path.join(path, prefix + '.stack.fits')
    #print(f_2d)
    twod_hdu = pyfits.open(f_2d)
    #twod_hdu = pyfits.open(os.path.join(path,f_2d))
    #print(twod_hdu.info())
    #print()

    if not fig:
        fig = plt.figure(figsize=(12,10))

    #           row column
    gs = GridSpec(3, 4, hspace=0.4)

    ax = {}

    ax[0] = fig.add_subplot(gs[:1, 0:2])
    ax[1] = fig.add_subplot(gs[:1, 2:])
    ax[2] = fig.add_subplot(gs[1, :])
    ax[3] = fig.add_subplot(gs[2, :])


    # Z_MAP, CHIMIN and DOF
    
    #ax[0].text(0.95, 0.96, root + '\n'+'ID={0:<5d}  z={1:.4f}'.format(id, zfit_head['z_map']), ha='right', va='top', transform=ax[0].transAxes, fontsize=9)
    ax[0].text(0.95, 0.96, prefix + '\n'+'z_fit={0:.4f}'.format(zfit_head['z_map']), ha='right', va='top', transform=ax[0].transAxes, fontsize=9)

    zmi, zma = zfit_stack['zgrid'].min(), zfit_stack['zgrid'].max()
    if (zma-zmi) > 5:
        ticks = np.arange(np.ceil(zmi), np.floor(zma)+0.5, 1)
        lz = np.log(1+zfit_stack['zgrid'])
        ax[0].plot(lz, np.log10(zfit_stack['pdf']), color='k')
        ax[0].set_xticks(np.log(1+ticks))
        ax[0].set_xticklabels(np.cast[int](ticks))
        ax[0].set_xlim(lz.min(), lz.max())
    else:
        ax[0].plot(zfit_stack['zgrid'], np.log10(zfit_stack['pdf']), color='k')
        ax[0].set_xlim(zmi, zma)

    ax[0].set_xlabel(r'$z$')
    ax[0].set_ylabel(r'$\log\ p(z)$'+' / '+ r'$\chi^2=\frac{{{0:.0f}}}{{{1:d}}}={2:.2f}$'.format(zfit_head['chimin'], zfit_head['DoF'], zfit_head['chimin']/zfit_head['DoF']))
    ax[0].set_yticks([1,4,9,16,25])

    pzmax = np.log10(zfit_stack['pdf'].max())
    ax[0].set_ylim(pzmax-6, pzmax+0.9)
    ax[0].grid()
    ax[0].yaxis.set_major_locator(MultipleLocator(base=1))
    
    
    
    #pz_max = np.log10(zfit_stack['pdf'].max())
    
    #ax[0].plot(zfit_stack['zgrid'], np.log10(zfit_stack['pdf']), label='Stacked')
    #ax[0].set_xlim(0.0, 3); 
    ##ax[0].semilogy(); 
    #ax[0].grid()
    ##ax[0].set_ylim(1.e-50, 1e4)
    #ax[0].set_ylim(pz_max-6, pz_max+0.9)
    #ax[0].set_xlabel('z'); ax[0].set_ylabel('PDF(z)'); #plt.legend()

    ax[1].plot(zfit_stack['zgrid'], zfit_stack['risk'], label='Stacked')
    ax[1].set_xlim(0.0, 3); ax[1].semilogy(); ax[1].grid()
    ax[1].set_xlabel('z'); ax[1].set_ylabel('risk'); #ax[2].legend()
    
    #print('Continuum template, cont1d: ', templ['continuum'].__class__)
    #print(templ.colnames)
    
    

    
    ax[2].plot(templ['wave']/1.e4, templ['continuum'], label='continuum')
    ax[2].plot(templ['wave']/1.e4, templ['full'], label='total')

    ax[2].grid()
    ax[2].set_xlabel(r'$\lambda$ (microns)')
    ax[2].set_ylabel(r'F$_\lambda$ (erg/s/cm$^2$/$\AA$)')
    #ax[2].legend()
    
    
        

    yy0 = 1
    yy1 = -1    

    for i,disp in enumerate(dispersers):

        flux = tbl[i]['flux']/tbl[i]['flat']
        wav = tbl[i]['wave']/1e4

        filt = (wav > w0) & (wav < w1)

        ax[3].errorbar(wav, flux, yerr=tbl[i]["err"]/tbl[i]['flat'], 
                    color=colors[i], marker='.', linestyle='None', alpha=0.5, label=disp.upper())
        
        if not i: label = "model"
        else: label = ""
        ax[3].plot(wav, tbl[i]['line']/tbl[i]['flat'], color="g", alpha=0.5, label=label)

        if not i: label = "contam"
        else: label = ""
        ax[3].plot(wav, tbl[i]['contam']/tbl[i]['flat'], color="r", alpha=0.5, label=label)

        if np.nanmin(flux[filt]) < yy0: yy0 = np.nanmin(flux[filt])
        if np.nanmax(flux[filt]) > yy1: yy1 = np.nanmax(flux[filt])
    

    #if y0 == -1: y0 = np.min(templ['continuum'] + templ['full'])
    #if y1 == -1: y1 = np.max(templ['continuum'] + templ['full'])
    if y0 == -1:
        if yy0 > 0:
            y0 = yy0/yscale
        else:
            y0 = yscale*yy0
    if y1 == -1: y1 = yscale*yy1


    if z_in:
        ax[3].text(0.95, 0.96, 'z_fit={0:.4f}'.format(zfit_head['z_map']), ha='right', va='top', transform=ax[3].transAxes, fontsize=9)
        ax[3].text(0.95, 0.90, 'z_input={0:.4f}'.format(z_in), color="b", ha='right', va='top', transform=ax[3].transAxes, fontsize=9)
        z0 = z_in
    else:
        ax[3].text(0.95, 0.96, 'z_fit={0:.4f}'.format(zfit_head['z_map']), color="b", ha='right', va='top', transform=ax[3].transAxes, fontsize=9)
        z0 = zfit_head['z_map']

    for line,wemit in emlines:
        wobs = (1+z0)*wemit
        #FWHM = vel * wobs / c_km
        if wobs/1e4 > w0 and wobs/1e4 < w1:

            ax[3].plot([wobs/1e4,wobs/1e4],[y0,y1],":",lw=lw,c="b")
        
            if labels:
                ax[3].text(wobs/1e4,0.7*y1,line,fontsize=fontsize,
                        rotation='vertical',
                        horizontalalignment='center',
                        verticalalignment='center')

    ax[2].xaxis.set_major_locator(MultipleLocator(0.1))
    ax[2].set_xlim(w0, w1); ax[2].set_ylim(y0,y1)


    ax[3].xaxis.set_major_locator(MultipleLocator(0.1))
    ax[3].set_xlim(w0, w1); ax[3].set_ylim(y0,y1)
    #ax[3].semilogy(); 
    ax[3].grid()
    ax[3].set_xlabel(r'$\lambda$ (microns)')
    ax[3].set_ylabel(r'F$_\lambda$ (erg/s/cm$^2$/$\AA$)')

    ax[3].legend()
    
    # Gabe's routine
    multifit.show_drizzle_HDU(twod_hdu, mask_segmentation=True, diff=True, average_only=False, cmap=cmap, scale_size=scale_size)


    return ax




emlines = [["OVI",         1038.0],         # 0
           ["Ly$\\alpha$", 1215.67],        # 1
           ["CIV",     1550.0],             # 2
           ["CIII]",   1909.],              # 3
           ["CII]",    2327.],              # 4
           ["MgII",    2796.4],             # 5
           ["MgII",    2803.5],             # 6
           ["NeV",     3326.],              # 7
           ["[OII]",   3727.],  # O2        # 8
           ["[NeIII]", 3868.7],             # 9
           ["H$\gamma$",  4340.5],  # Hg    # 10
           ["[OIII]",  4363.0],  # O31      # 11
           ["H$\\beta$",   4861.3],  # Hb   # 12
           ["[OIII]",  4959.0],  # O32      # 13
           ["[OIII]",  5007.0],  # O33      # 14
           ["[NII]",   6548.1],             # 15
           ["H$\\alpha$",  6562.8],  # Ha   # 16
           ["[NII]",   6583.0],             # 17
           ["[SII]",   6717.0],             # 18
           ["[SII]",   6731.0],             # 19
           ["[SIII]",  9069.0],
           ["[SIII]",  9532.0],
           ["P$\\delta$", 10049.8],  # Pd   # 20
           ["P$\\gamma$", 10938.0],  # Pg   # 21
           ["P$\\beta$",  12818.1],  # Pb   # 22
           ["P$\\alpha$", 18750.1],  # Pa   # 23 
           ["Br$\\delta$", 19440.0],  # Br-d (wikipedia, not exact)
           ["Br$\\gamma$", 21660.0],  # Br-g (wikipedia, not exact)
           ["Br$\\beta$",  26250.0],  # Br-b (wikipedia, not exact)
           ["Br$\\alpha$", 40510.0],  # Br-a (wikipedia, not exact) 
          ]

# http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/hydspec.html
# http://articles.adsabs.harvard.edu//full/1934ApJ....80...19M/0000022.000.html


def redshift_analysis_od(root,id,save_lines=3):

    full_obj = OrderedDict()
   
    full_hdu = pyfits.open('{0}_{1:05d}.full.fits'.format(root, id))
    h0 = full_hdu[0].header
    line_keys = ['ID','RA','DEC','REDSHIFT','NUMLINES']
    #print(h0)
    for lk in line_keys:
        #h0[lk],h0.comments[lk]
        #print(h0[lk],)
        full_obj[lk] = h0[lk]
    
    # read all of the lines
    #print(h0['NUMLINES'])
    lines_list = []
    for i in range(h0['NUMLINES']):
        #print(i+1)
        od = OrderedDict({"line":h0["LINE%03d" % (i+1)],
                          "flux":h0["FLUX%03d" % (i+1)],
                          "err":h0["ERR%03d" % (i+1)]})                  
        lines_list.append(od)
                          
                          
    #print(lines_list)
    line_tbl = Table(lines_list)
    line_tbl.sort("flux")
    line_tbl.reverse()
    
    # store only 3 (for now)
    for i in range(save_lines):
        try:
            full_obj["LINE%03d" % (i+1)] = line_tbl['line'][i]
            full_obj["FLUX%03d" % (i+1)] = line_tbl['flux'][i]
            full_obj["ERR%03d" % (i+1)] = line_tbl['err'][i]
        except IndexError:
            full_obj["LINE%03d" % (i+1)] = "None"
            full_obj["FLUX%03d" % (i+1)] = 99.9
            full_obj["ERR%03d" % (i+1)] = -99.9
    
    
    h1 = full_hdu['ZFIT_STACK'].header
    zfit_keys = ['CHIMIN','DOF','Z02','Z16','Z50','Z84','Z97','ZWIDTH1',
                 'ZWIDTH2','Z_RISK','MIN_RISK','Z_MAP','GAM_LOSS']
    
    for zfk in zfit_keys:
        #h1[zfk],h1.comments[zfk]
        #print(h1[zfk],)
        full_obj[zfk] = h1[zfk]
    
    return full_obj

def plot_redshifts(tbl):

    fig = plt.figure()
    p1 = fig.add_subplot(111)
    p1.scatter(tbl["ZWIDTH1"], tbl["ZWIDTH2"])

    p1.plot([0.0001,10.0],[0.05,0.05],"--",c="k")
    p1.plot([0.05,0.05],[0.0001,10.0],"--",c="k")

    p1.set_xlabel("zwidth1 [16th and 84th p(z) percentile]")
    p1.set_ylabel("zwidth2 [2.5th and 97.5th p(z) percentile]")
    p1.set_xscale("log")
    p1.set_yscale("log")
    p1.set_xlim(0.0001,10.0)
    p1.set_ylim(0.0001,10.0)

    fig = plt.figure(figsize=(12,4))
    p1 = fig.add_subplot(121)
    p1.scatter(tbl["REDSHIFT"], tbl["ZWIDTH1"])
    #p1.plot([0,3.5],[0.1,0.1],"--",c="k")
    p1.plot([0,3.5],[0.05,0.05],"--",c="k")
    p1.set_xlabel("Redshift")
    p1.set_ylabel("zwidth1 [16th and 84th p(z) percentile]")
    p1.set_yscale("log")
    p1.set_ylim(0.0001,10.0)

    p2 = fig.add_subplot(122)
    p2.scatter(tbl["REDSHIFT"], tbl["ZWIDTH2"])
    #p2.plot([0,3.5],[0.1,0.1],"--",c="k")
    p2.plot([0,3.5],[0.05,0.05],"--",c="k")
    p2.set_xlabel("Redshift")
    p2.set_ylabel("zwidth2 [2.5th and 97.5th p(z) percentile]")
    p2.set_yscale("log")
    p2.set_ylim(0.0001,10.0)
    
    print(len(tbl))

def plot_redshifts2(tbl,zcut=0.01):

    filt1 = tbl["ZWIDTH1"] < zcut
    tbl_filt1 = tbl[filt1]

    fig = plt.figure(figsize=(12,10))

    p1 = fig.add_subplot(221)
    p1.scatter(tbl["REDSHIFT"], tbl["ZWIDTH1"], s=2, c="k")
    p1.scatter(tbl_filt1["REDSHIFT"], tbl_filt1["ZWIDTH1"])
    #p1.plot([0,3.5],[0.1,0.1],"--",c="k")
    p1.plot([0,3.5],[zcut,zcut],"--",c="k")
    p1.set_xlabel("Redshift")
    p1.set_ylabel("z$_{width1}$ [16th and 84th p(z) percentile]")
    p1.set_yscale("log")
    p1.set_ylim(0.0001,10.0)

    p2 = fig.add_subplot(222)
    p2.scatter(tbl["REDSHIFT"], tbl["MAG_AUTO"], s=2, c="k")
    p2.scatter(tbl_filt1["REDSHIFT"], tbl_filt1["MAG_AUTO"])
    #p2.plot([0,3.5],[0.1,0.1],"--",c="k")
    #p2.plot([0,3.5],[0.05,0.05],"--",c="k")
    p2.set_xlabel("Redshift")
    p2.set_ylabel("Mag")
    #p2.set_yscale("log")
    #p2.set_ylim(0.0001,10.0)

    filt2 = tbl["ZWIDTH2"] < zcut
    tbl_filt2 = tbl[filt2]

    p3 = fig.add_subplot(223)
    p3.scatter(tbl["REDSHIFT"], tbl["ZWIDTH2"], s=2, c="k")
    p3.scatter(tbl_filt2["REDSHIFT"], tbl_filt2["ZWIDTH2"])
    #p3.plot([0,3.5],[0.1,0.1],"--",c="k")
    p3.plot([0,3.5],[zcut,zcut],"--",c="k")
    p3.set_xlabel("Redshift")
    p3.set_ylabel("z$_{width2}$ [2.5th and 97.5th p(z) percentile]")
    p3.set_yscale("log")
    p3.set_ylim(0.0001,10.0)

    p4 = fig.add_subplot(224)
    p4.scatter(tbl["REDSHIFT"], tbl["MAG_AUTO"], s=2, c="k")
    p4.scatter(tbl_filt2["REDSHIFT"], tbl_filt2["MAG_AUTO"])
    #p4.plot([0,3.5],[0.1,0.1],"--",c="k")
    #p4.plot([0,3.5],[0.05,0.05],"--",c="k")
    p4.set_xlabel("Redshift")
    p4.set_ylabel("Mag")
    #p4.set_yscale("log")
    #p4.set_ylim(0.0001,10.0)

    print(len(tbl))
    print(len(tbl[filt1]))
    print(len(tbl[filt2]))
    
    return filt1,filt2

###################################
###################################
###################################
# Function to fit a single redshift
###################################
###################################
###################################
    
#def grizli_fit_z(sim,root,id,fwhm,temp0,temp1):
def grizli_fit_z(root,id,fwhm,z1=0.05,z2=3.0):

    #collected = gc.collect()
    #print("Garbage collector: collected %d objects." % (collected))

    print(id)
    
    #beams = OrderedDict()
    
    is_cgs, spectrum_1d, b = sim.object_dispersers[id]
    cutout = grizli.model.BeamCutout(sim, b['A'], min_sens=0,) # min_mask=0) 
    
    cutout.beam.compute_model()  
    cutout.contam = cutout.beam.cutout_from_full_image(sim.model)
    if id in sim.object_dispersers:
        cutout.contam -= cutout.beam.model
    
    #beams[sim.grism.instrument] = cutout
    
    ###########################
    # new lines of code to test
    ###########################
    hdu = cutout.write_fits(get_hdu=True)
    ext = 0
    hdu[ext].header['EXPTIME'] = hdu['SCI'].header['EXPTIME']
    
    beam = 'beam_%05d.grism.A.fits' % (id)
    hdu.writeto(beam,clobber=True)
    
    mb = multifit.MultiBeam([beam], fcontam=0.2, group_name=root, psf=False, min_sens=0.05)
    mb.write_master_fits()
    
    # kludge
    os.remove(beam)
    ###########################
    
    fitting.run_all(id, temp0, temp1, fit_only_beams=True, fwhm=fwhm, zr=[z1, z2], 
                    dz=[0.004, 0.0002], fitter=['nnls', 'bounded'], group_name=root)

    collected = gc.collect()
    print("Garbage collector: collected %d objects." % (collected))

    return 1

####################################################
####################################################
####################################################
# Function to fit all redshifts calling grizli_fit_z
####################################################
####################################################
####################################################

def grizli_fit_all_z_3(HOME_PATH='/Users/gwalth/data/Roman/grizli/',
                       root='my_roman_sims',
                       field='field_0001',
                       pickle_file='Roman_GrismFLT.pickle',
                       sources = [],
                       fwhm=325, mag_limit=30,z1=0.05,z2=3.0):


    memory_tracker = t##racker.SummaryTracker()
    #memory_tracker.print_diff()

    print('HOME_PATH = ', HOME_PATH)

    global sim
    with open(pickle_file, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        sim = pickle.load(f)

    if not sources: 
        #phot = Table.read(phot_file,format="ascii")
        #sources = [id for id in phot['NUMBER']]

        all_src,magcut_src,extract_src = check_sims(sim, mag_limit)
        sources = [id for id in extract_src['NUMBER']]
    #print(sources)

    print("Computing redshifts for %i sources" % (len(sources)))

    os.chdir(os.path.join(HOME_PATH, root, 'Extraction',field))
    
    # Fitting templates
    global temp0
    global temp1
    # First is set with combined emission line complexes for the redshift fit 
    # (don't allow infinite freedom) of the line ratios / fluxes
    temp0 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=True, stars=False, 
                                         full_line_list=None,  continuum_list=None, 
                                         fsps_templates=True)
    
    # Second set has individual line templates for fitting the line fluxes
    temp1 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=False, stars=False, 
                                         full_line_list=None, continuum_list=None, 
                                         fsps_templates=True)
    


    t0 = time.time()

    for src_id in sources:
        #grizli_fit_z(sim,root,src_id,fwhm,temp0,temp1)
        grizli_fit_z(root,src_id,fwhm,z1=z1,z2=z2)
        #memory_tracker.print_diff()

        #print()
        #all_objects = muppy.get_objects()
        #print("All objects = %i" % (len(all_objects)))
        #sum1 = summary.summarize(all_objects)
        #summary.print_(sum1)
        memory_tracker.print_diff()

    t1 = time.time()
    print("%.2f seconds" % (t1-t0))



def grizli_fit_all_z_2(HOME_PATH='/Users/gwalth/data/Roman/grizli/',
                       root='my_roman_sims',
                       field='field_0001',
                       pickle_file='Roman_GrismFLT.pickle',
                       sources = [],
                       fwhm=325, mag_limit=30,z1=0.05,z2=6.0):

    print('HOME_PATH = ', HOME_PATH)

    global sim
    with open(pickle_file, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        sim = pickle.load(f)

    if not sources: 
        #phot = Table.read(phot_file,format="ascii")
        #sources = [id for id in phot['NUMBER']]

        all_src,magcut_src,extract_src = check_sims(sim, mag_limit)
        sources = [id for id in extract_src['NUMBER']]
    #print(sources)

    print("Computing redshifts for %i sources" % (len(sources)))

    os.chdir(os.path.join(HOME_PATH, root, 'Extraction',field))
    
    # Fitting templates
    global temp0
    global temp1
    # First is set with combined emission line complexes for the redshift fit 
    # (don't allow infinite freedom) of the line ratios / fluxes
    temp0 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=True, stars=False, 
                                         full_line_list=None,  continuum_list=None, 
                                         fsps_templates=True)
    
    # Second set has individual line templates for fitting the line fluxes
    temp1 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=False, stars=False, 
                                         full_line_list=None, continuum_list=None, 
                                         fsps_templates=True)
    


    t0 = time.time()

    for src_id in sources:
        #grizli_fit_z(sim,root,src_id,fwhm,temp0,temp1)
        grizli_fit_z(root,src_id,fwhm,z1=z1,z2=z2)

    t1 = time.time()
    print("%.2f seconds" % (t1-t0))


def grizli_fit_all_z(HOME_PATH='/Users/gwalth/data/Roman/grizli/',
                     root='my_roman_sims',
                     pickle_file='Roman_GrismFLT.pickle',
                     sources = [],
                     fwhm=325, mag_limit=30,z1=0.05,z2=3.0):

    print('HOME_PATH = ', HOME_PATH)

    global sim
    with open(pickle_file, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        sim = pickle.load(f)

    if not sources: 
        #phot = Table.read(phot_file,format="ascii")
        #sources = [id for id in phot['NUMBER']]

        all_src,magcut_src,extract_src = check_sims(sim, mag_limit)
        sources = [id for id in extract_src['NUMBER']]
    #print(sources)

    print("Computing redshifts for %i sources" % (len(sources)))

    os.chdir(os.path.join(HOME_PATH, root, 'Extraction'))
    
    # Fitting templates
    global temp0
    global temp1
    # First is set with combined emission line complexes for the redshift fit 
    # (don't allow infinite freedom) of the line ratios / fluxes
    temp0 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=True, stars=False, 
                                         full_line_list=None,  continuum_list=None, 
                                         fsps_templates=True)
    
    # Second set has individual line templates for fitting the line fluxes
    temp1 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=False, stars=False, 
                                         full_line_list=None, continuum_list=None, 
                                         fsps_templates=True)
    


    t0 = time.time()

    for src_id in sources:
        #grizli_fit_z(sim,root,src_id,fwhm,temp0,temp1)
        grizli_fit_z(root,src_id,fwhm,z1=z1,z2=z2)

    t1 = time.time()
    print("%.2f seconds" % (t1-t0))


def grizli_fit_all_z_async_2(HOME_PATH='/Users/gwalth/data/Roman/grizli/',
                             root='my_roman_sims',
                             pickle_file='Roman_GrismFLT.pickle',
                             field='field_0001',
                             sources = [],
                             N_CPU='',
                             fwhm=325, mag_limit=30,z1=0.05,z2=6.0):

    print('HOME_PATH = ', HOME_PATH)

    global sim
    with open(pickle_file, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        sim = pickle.load(f)

    if not sources:
        #phot = Table.read(phot_file,format="ascii")
        #sources = [id for id in phot['NUMBER']]

        all_src,magcut_src,extract_src = check_sims(sim, mag_limit)
        sources = [id for id in extract_src['NUMBER']]

    print("Computing redshifts for %i sources" % (len(sources)))

    #os.chdir(os.path.join(HOME_PATH, root, 'Extraction'))
    os.chdir(os.path.join(HOME_PATH, root, 'Extraction',field))

    # Fitting templates

    global temp0
    global temp1
    # First is set with combined emission line complexes for the redshift fit 
    # (don't allow infinite freedom) of the line ratios / fluxes
    temp0 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=True, stars=False,
                                        full_line_list=None,  continuum_list=None,
                                        fsps_templates=True)

    # Second set has individual line templates for fitting the line fluxes
    temp1 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=False, stars=False,
                                        full_line_list=None, continuum_list=None,
                                        fsps_templates=True)



    # version 2 
    t0 = time.time()
    # Step 1: Init multiprocessing.Pool()
    print()
    print("Number of CPUs =",mp.cpu_count())
    if not N_CPU:
        N_CPU = mp.cpu_count()
    print("Using %i CPUs" % N_CPU)
    #pool = mp.Pool(mp.cpu_count())
    pool = mp.Pool(N_CPU)

    processes = []
    # Step 3: Use loop to parallelize
    for src_id in sources:
        p = pool.apply_async(grizli_fit_z, args=(root,src_id,fwhm))
        processes.append(p)

    pool.close()
    pool.join()

    result = [p.get() for p in processes]

    t1 = time.time()
    print("%.2f seconds" % (t1-t0))
    



def grizli_fit_all_z_async(HOME_PATH='/Users/gwalth/data/Roman/grizli/',
                           root='my_roman_sims',
                           pickle_file='Roman_GrismFLT.pickle',
                           phot_file = '',
                           sources = [],
                           fwhm=325,
                           N_CPU=''):

    print('HOME_PATH = ', HOME_PATH)

    if not sources: 
        phot = Table.read(phot_file,format="ascii")
        sources = [id for id in phot['NUMBER']]

    print("Computing redshifts for %i sources" % (len(sources)))

    global sim
    with open(pickle_file, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        sim = pickle.load(f)



    os.chdir(os.path.join(HOME_PATH, root, 'Extraction'))
    
    # Fitting templates
    
    global temp0
    global temp1
    # First is set with combined emission line complexes for the redshift fit 
    # (don't allow infinite freedom) of the line ratios / fluxes
    temp0 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=True, stars=False, 
                                        full_line_list=None,  continuum_list=None, 
                                        fsps_templates=True)
    
    # Second set has individual line templates for fitting the line fluxes
    temp1 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=False, stars=False, 
                                        full_line_list=None, continuum_list=None, 
                                        fsps_templates=True)



    # version 2 
    t0 = time.time()
    # Step 1: Init multiprocessing.Pool()
    if not N_CPU:
        N_CPU = mp.cpu_count()
    print("Using %i CPUs" % N_CPU)
    #pool = mp.Pool(mp.cpu_count())
    pool = mp.Pool(N_CPU)

    results = []
    # Step 3: Use loop to parallelize
    for src_id in sources:
        #pool.apply_async(grizli_fit_z, args=(sim,root,src_id,fwhm,temp0,temp1),
        #                 callback=collect_result)
        pool.apply_async(grizli_fit_z, args=(root,src_id,fwhm),
                         callback=collect_result)

    # Step 4: Close Pool and let all the processes complete    
    pool.close()
    pool.join()  # postpones the execution of next line of code until all processes in the queue are done.

    t1 = time.time()
    print("%.2f seconds" % (t1-t0))


def grizli_fit_all_z_sync(HOME_PATH='/Users/gwalth/data/Roman/grizli/',
                          root='my_roman_sims',
                          pickle_file='Roman_GrismFLT.pickle',
                          phot_file = '',
                          sources = [],
                          fwhm=325,
                          N_CPU=''):



    print('HOME_PATH = ', HOME_PATH)

    if not sources: 
        phot = Table.read(phot_file,format="ascii")
        sources = [id for id in phot['NUMBER']]
    
    print("Computing redshifts for %i sources" % (len(sources)))

    global sim
    with open(pickle_file, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        sim = pickle.load(f)



    os.chdir(os.path.join(HOME_PATH, root, 'Extraction'))
    
    # Fitting templates
    
    global temp0
    global temp1
    # First is set with combined emission line complexes for the redshift fit 
    # (don't allow infinite freedom) of the line ratios / fluxes
    temp0 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=True, stars=False, 
                                        full_line_list=None,  continuum_list=None, 
                                        fsps_templates=True)
    
    # Second set has individual line templates for fitting the line fluxes
    temp1 = grizli.utils.load_templates(fwhm=fwhm, line_complexes=False, stars=False, 
                                        full_line_list=None, continuum_list=None, 
                                        fsps_templates=True)
    


    # version 1 
    t0 = time.time()
    # Step 1: Init multiprocessing.Pool()
    if not N_CPU:
        N_CPU = mp.cpu_count()
    print("Using %i CPUs" % N_CPU)
    #pool = mp.Pool(mp.cpu_count())
    pool = mp.Pool(N_CPU)

    # Step 2: `pool.apply` the `grizli_fit_z`
    #[pool.apply(grizli_fit_z, args=(sim,root,src_id,fwhm,temp0,temp1)) for src_id in sources]
    [pool.apply(grizli_fit_z, args=(root,src_id,fwhm)) for src_id in sources]

    # Step 3: Don't forget to close
    pool.close()
    t1 = time.time()
    print("%.2f seconds" % (t1-t0))

########################
########################
########################
# Forking test functions
########################
########################
########################

def howmany_within_range(row, minimum, maximum):
    """Returns how many numbers lie within `maximum` and `minimum` in a given `row`"""
    count = 0
    for n in row:
        if minimum <= n <= maximum:
            count = count + 1
    return count

def test_func():

    # Prepare data
    np.random.RandomState(100)
    arr = np.random.randint(0, 10, size=[200000, 5])
    data = arr.tolist()
    data[:5]

    t0 = time.time()
    results = [howmany_within_range(row, 4, 8) for row in data]
    t1 = time.time()
    print(t1-t0,"seconds")
    print(results[:10])

def test_func_sync():
    # Prepare data
    np.random.RandomState(100)
    arr = np.random.randint(0, 10, size=[200000, 5])
    data = arr.tolist()
    data[:5]

    t0 = time.time()
    # Step 1: Init multiprocessing.Pool()
    pool = mp.Pool(mp.cpu_count())

    # Step 2: `pool.apply` the `howmany_within_range()`
    results = [pool.apply(howmany_within_range, args=(row, 4, 8)) for row in data]

    # Step 3: Don't forget to close
    pool.close()    

    t1 = time.time()
    print(t1-t0,"seconds")
    print(results[:10])

def collect_result(result):
    global results
    results.append(result)

########################
# other useful functions
########################


def get_sim_info(pickle_file='Roman_GrismFLT.pickle', mag_limit=30):


    with open(pickle_file, 'rb') as f:
        # The protocol version used is detected automatically, so we do not
        # have to specify it.
        sim = pickle.load(f)

    all_src,magcut_src,extract_src = check_sims(sim, mag_limit)
    sources = [id for id in extract_src['NUMBER']]

    return sources
#print(Roman)
#

## ~ 5GB file
#with open('Roman_GrismFLT.pickle', 'wb') as f:
##    # Pickle the 'data' dictionary using the highest protocol available.
#    pickle.dump(Roman, f, pickle.HIGHEST_PROTOCOL)
##
#
#
#print(Roman)
#print(Roman.__dict__)

#/local/RomanSims/grizli/sims/sim_v3/Prep/field_0002
#
#default.nnw                                             
#Euclid_Roman_4deg2_field_0002_v3_direct.fits
#Euclid_Roman_4deg2_field_0002_v3_slitless.fits
#Euclid_Roman_4deg2_field_0002.fits
#
#Roman_GrismFLT.pickle
#Roman.param
#Roman.sex
#
#/local/RomanSims/grizli/sims/sim_v3/Extraction
#
#/home/gwalth/data/Roman/grizli/sims/sim_v3


def prep_setup(input_path = "/home/gwalth/data/Roman/Galacticus/fields/",
               axe_path = "/home/gwalth/aXeSIM_Roman/OUTSIM/",
               PREP_PATH = "/home/gwalth/data/Roman/grizli/sims/",
               EXTRACTION_PATH = "/local/RomanSims/grizli/sims/",
               root = "sim_v3",
               fmt = "field_%04d",
               ver = "v3",
               field_start=20, field_end=23, field_list=[]):

    # SExtractor only works on cygnusc

    #input_path = "/home/gwalth/data/Roman/Galacticus/fields/"
    #input_path = "/local/RomanSims/Galacticus/fields/"
    #input_path = "/home/gwalth/data/Roman/grizli/sims/sim_v3/"
    #axe_path = "/home/gwalth/aXeSIM_Roman/OUTSIM/"
    #PREP_PATH = "/home/gwalth/data/Roman/grizli/sims/"
    #PREP_PATH = "/local/RomanSims/grizli/sims/"
    #EXTRACTION_PATH = "/home/gwalth/data/Roman/grizli/sims/"
    #EXTRACTION_PATH = "/local/RomanSims/grizli/sims/"
    #root = "sim_v3"
    #ver = "v3"

    sim_path = os.path.join(PREP_PATH, root, 'Prep')

    if not field_list:
        field_list = range(field_start,field_end+1,1)

    for i in field_list:

        field = fmt % (i)
        print("Preparing %s" % (field))
        print()

        cat = "Euclid_Roman_4deg2_%s.fits" % (field)
        spec = "Euclid_Roman_4deg2_%s_%s_slitless.fits" % (field, ver)
        img = "Euclid_Roman_4deg2_%s_%s_direct.fits" % (field, ver)

  
        prep_field_path = os.path.join(PREP_PATH, root, 'Prep', field)

        if not os.path.exists(prep_field_path):
            os.mkdir(prep_field_path)
            print("Creating directory:\n %s" % (prep_field_path))
            print()

        #os.chdir(prep_field_path)

        print("Creating symbolic links to:")
        if os.path.islink(prep_field_path + "/" + cat):
            os.remove(prep_field_path + "/" + cat)
            print("Removing old symbolic link")
        os.symlink(input_path + cat, prep_field_path + "/" + cat)
        print(cat)

        if os.path.islink(prep_field_path + "/" + spec):
            os.remove(prep_field_path + "/" + spec)
            print("Removing old symbolic link")
        os.symlink(axe_path + spec, prep_field_path + "/" + spec)
        print(spec)

        if os.path.islink(prep_field_path + "/" + img):
            os.remove(prep_field_path + "/" + img)
            print("Removing old symbolic link")
        os.symlink(axe_path + img, prep_field_path + "/" + img)
        print(img)
        print()

        print("Copying SExtractor parameter files:")
        shutil.copyfile(sim_path + "/Roman.param", prep_field_path + "/Roman.param")
        print("Roman.param")
        shutil.copyfile(sim_path + "/Roman.sex", prep_field_path + "/Roman.sex")
        print("Roman.sex")
        shutil.copyfile(sim_path + "/default.nnw", prep_field_path + "/default.nnw")
        print("default.nnw")
        print()
        print()

        #extract_field_path = os.path.join(EXTRACTION_PATH, root, 'Extraction', field)
        #if not os.path.exists(extract_field_path):
        #    os.mkdir(extract_field_path)
        #    print("Creating directory %s" % (extract_field_path))
        #    print()


#def extract_setup():

def add_noise(hdu, scale=1.0, background=0.5, exptime=1.e4, nexp=10, readnoise=10., img_ext='SCI', seed=None):

    head = hdu[img_ext].header
    naxis = (head['NAXIS1'], head['NAXIS2'])

    # Simple error model of read noise and sky background
    var = nexp*readnoise**2 + background*exptime

    # electrons / s
    rms = np.sqrt(var)/exptime

    poisson = hdu['ERR'].data*exptime

    hdu['ERR'].data = np.sqrt(poisson**2 + var)/exptime


    for name in [img_ext, 'ERR', 'DQ']:
        hdu[name].header['EXPTIME'] = exptime
        hdu[name].header['NEXP'] = nexp
        hdu[name].header['BUNIT'] = 'ELECTRONS/S'
        hdu[name].header['BACKGR'] = background
        hdu[name].header['CALCRMS'] = rms, 'Variance used for random noise'

    if seed is not None:
        np.random.seed(seed)
        hdu['ERR'].header['SEED'] = seed, 'Random seed'

    hdu[img_ext].data *= scale
    hdu[img_ext].data += np.random.normal(size=np.array(naxis).T)*rms

    return hdu

def add_noise2(hdu, background=0.5, exptime=1.e4, readnoise=10., seed=None):
    # background [e-/s]
    # readnoise [e-]

    head = hdu['SCI'].header
    naxis = (head['NAXIS1'], head['NAXIS2'])
    
    # N_src  # count rate from src [e-/s]
    # N_sky  # count rate from sky [e-/s]
    # N_dark # dark current
    # N_RN   # readnoise
    # t      # exposure time
    # npix   # number of pixels in aperture
    
    # signal = N_src*t
    # noise = np.sqrt(N_src*t + npix*(N_sky*t + N_dark*t + N_RN**2))
    
    for name in ['SCI', 'ERR', 'DQ']:
        hdu[name].header['EXPTIME'] = exptime
        #hdu[name].header['NEXP'] = nexp
        hdu[name].header['BUNIT'] = 'ELECTRONS/S'
        hdu[name].header['BACKGR'] = background
        #hdu[name].header['CALCRMS'] = rms, 'Variance used for random noise'

    if seed is not None:
        np.random.seed(seed)
        hdu['ERR'].header['SEED'] = seed, 'Random seed'
    
    signal = hdu['SCI'].data                                    # [e-]
    noise = np.sqrt(signal*exptime + background*exptime + readnoise**2) # [e-]
    #noise = signal + background*exptime + readnoise**2 # [e-]
    #noise = np.sqrt(background*exptime + readnoise**2)         # [e-]
    
    hdu['ERR'].data = noise/exptime # [e-/s]
    
    total_signal = np.random.poisson(lam=noise, size=noise.shape).astype("float64") # [e-]
    #total_signal = signal + noise
    final_signal = total_signal/exptime # [e-/s]

    hdu['SCI'].data = final_signal

    
    return hdu

def add_noise3(hdu, background=0.5, exptime=1.e4, readnoise=10., seed=None):
    # background [e-/s]
    # readnoise [e-]

    head = hdu['SCI'].header
    naxis = (head['NAXIS1'], head['NAXIS2'])
    
    # N_src  # count rate from src [e-/s]
    # N_sky  # count rate from sky [e-/s]
    # N_dark # dark current
    # N_RN   # readnoise
    # t      # exposure time
    # npix   # number of pixels in aperture
    
    # signal = N_src*t
    # noise = np.sqrt(N_src*t + npix*(N_sky*t + N_dark*t + N_RN**2))
    
    for name in ['SCI', 'ERR', 'DQ']:
        hdu[name].header['EXPTIME'] = exptime
        #hdu[name].header['NEXP'] = nexp
        hdu[name].header['BUNIT'] = 'ELECTRONS/S'
        hdu[name].header['BACKGR'] = background
        #hdu[name].header['CALCRMS'] = rms, 'Variance used for random noise'

    if seed is not None:
        np.random.seed(seed)
        hdu['ERR'].header['SEED'] = seed, 'Random seed'
    
    signal = hdu['SCI'].data                                    # [e-]
    noise = np.sqrt(signal*exptime + background*exptime + readnoise**2) # [e-]
    #noise = signal + background*exptime + readnoise**2 # [e-]
    #noise = np.sqrt(background*exptime + readnoise**2)         # [e-]

    # Simple error model of read noise and sky background
    var = nexp*readnoise**2 + background*exptime
    
    hdu['ERR'].data = noise/exptime # [e-/s]
    
    total_signal = np.random.poisson(lam=noise, size=noise.shape).astype("float64") # [e-]
    #total_signal = signal + noise
    final_signal = total_signal/exptime # [e-/s]

    hdu['SCI'].data = final_signal

    
    return hdu


def fake_noise(hdu, rms, exptime=1.e4):
    
    head = hdu['SCI'].header
    naxis = (head['NAXIS1'], head['NAXIS2'])
    
    signal = hdu['SCI'].data # [e-] very close to the counts
    
    noise = np.random.normal(size=np.array(naxis).T)*rms # [e-/s]

    hdu['ERR'].data = noise
    hdu['SCI'].data = signal/exptime+noise
    
    return hdu



def fake_euclid_direct(det_dict, flux_key="TU_FNU_H_NISP_MAG", mag_key="NIR_H", gain=1.0, verb=1,
                       collect_area = 9926., wav_cen = 17714., wav_width = 4999., background=1.12, 
                       exptime=100., nexp=1, readnoise=6., eff_tot = 0.8, heads=None):

    # VIS  0.101 "/pixel  -  6x6 CCDs (4kx4k pixels each) - 12 um
    # NISP 0.3   "/pixel  -  4x4 HgCdTe NIR detectors (2kx2k pixels each) - 18 um
    
    # det_tbl expected input is in Jy
    # output of direct image should be e-/s
    
    # erg/s/cm^s/Hz * nu * area / (h*nu) = erg/s/cm^s/Hz * area / h = photons/s/cm^2
    #flux_conv = Jy/h/gain # converts Jy to e-/s/cm^2
    #Apix = (pixel_size*mu)**2 # area of pixel [cm^2]
    #print(flux_conv)
    #print(Apix)
    #print(collect_area)
        
    all_direct = []
    
    fluxes = []    
    mags = []
    
    thumb_temp = 'NIS_catalog_file_%s.thm.beamA.fits'

    for i in range(0,4): 
        for j in range(0,4):
            print(i+1,j+1)
            num = '%i%i' % (i+1,j+1)

            head1 = heads['DET%s.SCI' % (num)]

            #gain = head1['GAIN_DET']
                
            direct = "Euclid_DET%s_direct.fits" % (num)
        
            all_direct.append(direct)

            data = np.zeros((2048,2048)) 
            err = np.zeros((2048,2048)) 
            dq = np.ones((2048,2048)) 

            #data += 1e-15

            hdus = pyfits.open(thumb_temp % (num))
            print(thumb_temp % (num))
            print(len(hdus))
            #print(hdu[0].header)
            #print(hdu['1191'].header)
            #print(hdu.info())
        
            det_tbl = det_dict[num]

            for k,hdu in enumerate(hdus[1:]):
                head = hdu.header
                thumb_norm = hdu.data
                
                yl,xl = thumb_norm.shape
                #box_hw_x = int(xl/2)
                #box_hw_y = int(yl/2)

                xpix,ypix = det_tbl['X_PIX','Y_PIX'][k]
            
                if xpix < 0 or xpix > 2048 or ypix < 0 or ypix > 2048: continue
    
                #xpix = int(xpix+0.5)
                #ypix = int(ypix+0.5)
            
                #y0 = ypix-box_hw_y
                #y1 = ypix+box_hw_y
                #x0 = xpix-box_hw_x
                #x1 = xpix+box_hw_x

                y0 = int(ypix - yl/2 + 0.5)
                y1 = int(ypix + yl/2 + 0.5)
                x0 = int(xpix - xl/2 + 0.5)
                x1 = int(xpix + xl/2 + 0.5)
                
                fnu_Jy = det_tbl[flux_key][k] # Jy                
                fnu = fnu_Jy*Jy # erg/s/cm^2/Hz
                flambda = fnu*Ang/((wav_cen*Ang)**2/c) # erg/s/cm^2/Ang
                E_phot = (h*c)/(wav_cen*Ang) # erg
                flux = flambda*wav_width*collect_area*eff_tot/E_phot # photons/s
                
                #######################
                            
                #Aobj = yl*xl*Apix
                #Aobj = collect_area
                #flux = det_tbl[flux_key][k]*flux_conv*Aobj*exptime*nexp
                fluxes.append(flux)
                mags.append(det_tbl[mag_key][k])
                flux_thumb = thumb_norm*flux
            
                if verb:
                    print(k)
                    print(head['EXTNAME'])
                    print(thumb_norm.shape)
                    print(flux)
                    print(np.sum(flux_thumb))
                    print(x0,x1)
                    print(y0,y1)
                    print()

        
                data[y0:y1,x0:x1] = data[y0:y1,x0:x1] + flux_thumb
                # e-


            hdu = pyfits.PrimaryHDU(header=head1)  
            hdu.name = 'PRIMARY'
    
            hdu_sci = pyfits.ImageHDU(data, header=head1)  
            hdu_sci.name = 'SCI'
            hdu_err = pyfits.ImageHDU(err, header=head1)
            hdu_err.name = 'ERR'
            hdu_dq = pyfits.ImageHDU(dq, header=head1)
            hdu_dq.name = 'DQ'

            hdu1 = pyfits.HDUList([hdu,hdu_sci,hdu_err,hdu_dq])
            hdu1.info()

            # hdu e-  --> e-/s
            #hdu2 = add_noise(hdu1, background=background, exptime=dir_exptime, 
            #                nexp=nexp, readnoise=readnoise)
            hdu2 = add_noise(hdu1, background=background, exptime=exptime, 
                              nexp=nexp, readnoise=readnoise)
            hdu2.writeto("../" + direct, overwrite=True)
            
            
    if verb:
        fig = plt.figure(figsize=(5,5))
        ax = fig.add_subplot(111)
        ax.hist(np.log10(fluxes),bins=30)
        ax.set_yscale("log")
        plt.show()
            
    return all_direct, fluxes, mags

############################################################################

def build_roman_model(HOME_PATH = '/home/gwalth/data/Roman/grizli/sims/',
                      root = "sim_v3", ver = "v3", verb=0,
                      mag_zero = 26.285, #mag_zero = 25.91, #SExtractor ZP
                      scale=1.0, background=1.12, nexp=1, readnoise=16,
                      direct_exptime=141.*8, slitless_exptime=301.*8,
                      field_start = 20, field_end = 23, field_list = [],
                      fmt = "field_%04d"):


    import time
    T0 = time.time()

    #import glob
    import os, sys
    #from collections import OrderedDict

    import matplotlib as mpl    
    import matplotlib.pyplot as plt
    #from matplotlib.gridspec import GridSpec
    #from matplotlib.ticker import MultipleLocator

    #from IPython.display import Image

    mpl.rcParams['figure.figsize'] = (10.0, 6.0)
    mpl.rcParams['font.size'] = 14
    mpl.rcParams['savefig.dpi'] = 72

    import numpy as np
    #from math import cos, sin, atan2, pi
    #from math import sqrt, log

    import astropy
    import astropy.io.fits as pyfits
    #import astropy.wcs as pywcs
    from astropy.table import Table #, join
    #from astropy.modeling import models

    #import drizzlepac
    #import photutils

    import grizli
    import grizli.model
    #import grizli.multifit
    #from grizli import utils, multifit, fitting
    #import grizli.fake_image
    #from grizli.pipeline import auto_script
    #from grizli import prep


    print('\n Python version: ', sys.version)
    print('\n Grizli version: ', grizli.__version__)
    print('\n Astropy version: ', astropy.__version__)
    print()



    print('HOME_PATH = ', HOME_PATH)
    print()
    #os.chdir(HOME_PATH)


    # loop over all fields

    #fields = range(field_start,field_end+1,1)
    #N = len(fields)
    #print("Modelling objects for %i fields" % (N))
    #print()


    if not field_list:
        field_list = range(field_start,field_end+1,1)

    N = len(field_list)
    print("Modelling objects for %i fields" % (N))
    print()

    for i in field_list:


    #for i in fields:

        field = fmt % (i)
        print("Building Roman model for %s" % (field))
        print()

        os.chdir(os.path.join(HOME_PATH, root, 'Prep',field))

        primer_cat = "Euclid_Roman_4deg2_%s.fits" % (field)
        slitless = "Euclid_Roman_4deg2_%s_%s_slitless.fits" % (field, ver)
        direct = "Euclid_Roman_4deg2_%s_%s_direct.fits" % (field, ver)

        print(direct)
        print(slitless)
    
        suffix = "_final"
        # updated FITS files
        direct_noise = direct.replace(".fits",suffix+".fits")
        slitless_noise = slitless.replace(".fits",suffix+".fits")
    
        # add noise to direct and slitless 
    
        ##############
        ### DIRECT ###
        ##############
        hdu = pyfits.open(direct)
        hdu.info()
        #del hdu[ext].header['']
    
        ext = 0
        hdu[ext].header['FILTER'] = 'H158'
        #hdu[ext].header['INSTRUME'] = 'WFIv1-GLW' # Pandeia's sensitivity function
        # conf version that works with both sims
        hdu[ext].header['INSTRUME'] = 'WFIv2-GLW' # Anahita's sensitivity function
        ##hdu[ext].header['EXTVER'] = ext
    
        ext = 1
        hdu[ext].header['FILTER'] = 'H158'
        #hdu[ext].header['INSTRUME'] = 'WFIv1-GLW' # Pandeia's sensitivity function
        hdu[ext].header['INSTRUME'] = 'WFIv2-GLW' # Anahita's sensitivity function
        hdu[ext].header['EXTVER'] = ext
    
        # Adding my noise
        #add_noise(hdu, scale=1.0, background=1.12, exptime=141., nexp=1, readnoise=16) # single exposure
        #add_noise(hdu, scale=1.0, background=1.12, exptime=141.*8, nexp=1, readnoise=16) # full depth 
        add_noise(hdu, scale=scale, background=background, exptime=direct_exptime,
                  nexp=nexp, readnoise=readnoise)
    
        #add_noise2(hdu, background=1.12, exptime=141., readnoise=16) # single exposure
        #add_noise2(hdu, background=1.12, exptime=141.*8, readnoise=16) # full depth
    
        #fake_noise(hdu, rms=0.01, exptime=141.)
        #fake_noise(hdu, rms=0.001, exptime=141.)
        hdu.writeto(direct_noise, overwrite=True, output_verify='fix')
        print("Writing",direct_noise)
    
    
    
        ###############
        ### SPECTRA ###
        ###############
        hdu = pyfits.open(slitless)
        hdu.info()
    
        ext = 0
        ##hdu[ext].header['FILTER'] = 'H158'
        #hdu[ext].header['INSTRUME'] = 'WFIv1-GLW' # Pandeia's sensitivity function
        hdu[ext].header['INSTRUME'] = 'WFIv2-GLW' # Anahita's sensitivity function
        hdu[ext].header['FILTER'] = 'GRISM'
        ##hdu[ext].header['EXTVER'] = ext
    
        ext = 1
        ##hdu[ext].header['FILTER'] = 'H158'
        #hdu[ext].header['INSTRUME'] = 'WFIv1-GLW' # Pandeia's sensitivity function
        hdu[ext].header['INSTRUME'] = 'WFIv2-GLW' # Anahita's sensitivity function
        hdu[ext].header['FILTER'] = 'GRISM'
        hdu[ext].header['EXTVER'] = ext
    
        # Adding my noise
        #add_noise(hdu, scale=1.0, background=1.12, exptime=301., nexp=1, readnoise=16) # single exposure
        #add_noise(hdu, scale=1.0, background=1.12, exptime=301.*8, nexp=1, readnoise=16) # full depth
        add_noise(hdu, scale=scale, background=background, exptime=slitless_exptime,
                  nexp=nexp, readnoise=readnoise)
    
        #add_noise2(hdu, background=1.12, exptime=301., readnoise=16) # single exposure
        #add_noise2(hdu, background=1.12, exptime=301.*8, readnoise=16) # full depth
    
        #fake_noise(hdu, rms=0.01, exptime=301.) # good test to detect Ha source at 
        #fake_noise(hdu, rms=0.001, exptime=301.) 
        hdu.writeto(slitless_noise, overwrite=True, output_verify='fix')
        print("Writing",slitless_noise)
    
    
    
        # Load primer for inspecting data quality
        primer_cat = "Euclid_Roman_4deg2_%s.fits" % (field)
        print(primer_cat)
        #primer = Table.read(primer_cat, format='ascii.sextractor')
        primer = Table.read(primer_cat, format='fits')
        primer[:10].show_in_notebook()
        print(len(primer))
        print(primer.colnames)
        
        primer.rename_column('num', 'NUMBER')
        primer.rename_column('x_pix', 'X_IMAGE')
        primer.rename_column('y_pix', 'Y_IMAGE')
        primer.rename_column('m_new', 'MAG_F1600W')
        
        
        
        if verb: 
            ### Show them!
            
            ncols = 2
            
            all = [direct,direct_noise,slitless,slitless_noise]
            print(all)
            
            #N = len(all)
            # not sure why this isn't working
            #nrows = -(-N/ncols) # returns ceiling of division
            #print(nrows)
            nrows = 2
            
            fig = plt.figure()
            
            for i,a in enumerate(all):
                pf = pyfits.open(a)
                image = pf['SCI'].data
                print(np.mean(image),np.median(image),np.std(image))
                
            
                ax = fig.add_subplot(nrows,ncols,i+1)
                ax.imshow(image, interpolation='Nearest', 
                          origin='lower', cmap='inferno',
                          vmin=-0.1,vmax=0.1)
                #ax.scatter(primer['X_IMAGE'], primer['Y_IMAGE'], s=100,
                #               edgecolor='green', facecolor='none', alpha=0.7)
            
            fig.tight_layout(pad=0.5)
            
            
            
            
            fig = plt.figure(figsize=(10,10))
            
            pf = pyfits.open(direct)
            image = pf['SCI'].data
            print(np.mean(image),np.median(image),np.std(image))
            
            #filt = primer['MAG_F1600W'] < 20
            filt = (primer['MAG_F1600W'] > 20.9) & (primer['MAG_F1600W'] < 21)
            
                
            
            ax = fig.add_subplot(111)
            ax.imshow(image, interpolation='Nearest', 
                       origin='lower', vmin=-0.1, vmax=0.5, cmap='inferno')
            ax.scatter(primer['X_IMAGE'][filt], primer['Y_IMAGE'][filt], s=100,
                        edgecolor='green', facecolor='none', alpha=0.7)
            for i in range(len(primer[filt])):
                ax.text(primer[filt]['X_IMAGE'][i], primer[filt]['Y_IMAGE'][i],"%.2f" % (primer[filt]['MAG_F1600W'][i]),c="green")
            
            fig.tight_layout(pad=0.5)
        
        
        
        # remove segmentation FITS to redo "next" step
        #!rm *_seg.fits
        #os.remove("*_seg.fits")
        
        print(direct_noise)
        
        ## Make SExtractor catalog
        prefix = direct_noise.replace(".fits","")
        
        sex = "Roman.sex"
                 
        cat = prefix + ".cat"   
        seg = prefix + "_seg.fits"
        bkg = prefix + "_bkg.fits"
        #aper = prefix + "_aper.fits"
        
        if not os.path.exists(seg):
            #os.system('wget http://www.stsci.edu/~brammer/grism/grizli_xdf_sextractor.tar.gz')
            #os.system('tar xzvf grizli_xdf_sextractor.tar.gz')
            
            #sex = prefix + ".sex"
            
            
            wht = direct_noise + "[2]"
            direct_ext = direct_noise + "[1]"
            
            checkimage_name = seg + "," + bkg
            #checkimage_name = seg + "," + bkg + "," + aper
            
            sex_str = 'sex ' + direct_ext + ' -c ' + sex + ' -WEIGHT_IMAGE ' + wht + \
                      ' -CHECKIMAGE_NAME ' + checkimage_name + ' -CATALOG_NAME ' + cat + \
                      ' -MAG_ZEROPOINT %.2f' % (mag_zero)
            os.system(sex_str)
            print(sex_str)
        
        # awk '{ printf "circle(%f, %f, 0.00007) # text={%.3f}\n", $4, $5, $42 }' GRS_FOV1_roll0_dx0_dy0_SCA1_direct.cat > GRS_FOV1_roll0_dx0_dy0_SCA1_direct.reg
        
        print(cat)
        

        # check pixel size of direct and slitless        
        wcs_pixel_scale(direct_noise)
        wcs_pixel_scale(slitless_noise)
        
        
        # Read SExtractor catalog
        print(cat)
        phot = Table.read(cat, format='ascii.sextractor') # ref_cat in multimission
        print(len(phot))
        print(np.min(phot["MAG_AUTO"]),np.max(phot["MAG_AUTO"]))


        if verb:
            fig = plt.figure()
            p = fig.add_subplot(111)
            p.hist(phot["MAG_AUTO"],range=(10,32),bins=44)
            p.set_xlabel("H158 Magnitude")
            p.set_ylabel("N")
            p.set_xlim(10,32)
        
        
        
        t0 = time.time()
        
        ### Roman GRS grism
        # allow simulation of objects at the edges
        #pad=0 # pixels
        pad = 800 # I think this may be optimal given the spectra size (only seems to add 10-20 spectra)
        
        # using the ZP Anihita gave me for the direct image
        
        # sims              old  new 
        #mag_limit = 18 #    34   11
        #mag_limit = 20 #   273   49
        #mag_limit = 22 #  1599  334 
        #mag_limit = 24 #  5177 1830
        #mag_limit = 26 # 10015 5704
        #mag_limit = 28 #       7172
        
        
        mag_limit = 30 
        
        #h, wcs = grizli.fake_image.wfirst_header(ra=ra, dec=dec, pa_aper=pa_aper, naxis=(2048,2048))
        #grizli.fake_image.make_fake_image(h, output='wfirst.fits', exptime=EXPTIME, nexp=NEXP)
        
        Roman = grizli.model.GrismFLT(grism_file=slitless_noise, verbose=True, pad=pad,  
                                      ref_file=direct_noise, ref_ext=1,
                                      seg_file=seg, shrink_segimage=True)
        
        Roman_cat = Roman.blot_catalog(phot, sextractor=True) 
        Roman.catalog = Roman_cat
        
        mask = Roman_cat['MAG_AUTO'] < mag_limit
        print('N=%d' %(mask.sum()))
        #Roman.compute_full_model(compute_beams=['A'], mask=mask, verbose=False)
        Roman.compute_full_model(ids=Roman_cat['NUMBER'][mask], mags=Roman_cat['MAG_AUTO'][mask])
        #Roman.compute_full_model(ids=Roman_cat['NUMBER'][mask], mags=Roman_cat['MAG_AUTO'][mask], verbose=True)
        
        t1 = time.time()
        print()
        print("Finished computing Roman model in %.1f seconds" % (t1-t0))
        print()
        
        
        
        
        t0 = time.time()
        print(Roman)
        # ~ 5GB file
        with open('Roman_GrismFLT.pickle', 'wb') as f:
            # Pickle the 'data' dictionary using the highest protocol available.
            pickle.dump(Roman, f, pickle.HIGHEST_PROTOCOL)
        t1 = time.time()
        
        print()
        print("Finished saving pickle in %.1f seconds" % (t1-t0))
        print()
        
        
        # check Roman objects    
        print(Roman)
        print(Roman.__dict__)

    print()        
    T1 = time.time()
    print("Total time to compute Roman models = %.1f seconds" % (T1-T0))
    print("Number of fields = %i" % (N))
    print()


##############################
##############################
##############################
# JWST NIRISS redshift fitting
##############################
##############################
##############################


def grizli_fit_all_z_niriss(HOME_PATH='/Volumes/data3/grizli/',
                       root='my_data',
                       sources = [],
                       fwhm=325, mag_limit=30,z1=0.05,z2=12.0):

    #os.chdir(os.path.join(HOME_PATH + root, 'Extractions'))

    global grp
    grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), 
                            catalog='{0}-ir.cat.fits'.format(root), 
                            cpu_count=-1, sci_extn=1, pad=800)
    print("Number of FLTs:")
    print(len(grp.FLTs))


    if not sources: 
        all_ids = [id for flt in grp.FLTs for id in flt.object_dispersers]
        sources = list(set(all_ids))
        sources.sort()

    print("Computing redshifts for %i sources" % (len(sources)))


    pline={'kernel': 'square', 'pixfrac': 0.5, 'pixscale': 0.04, 'size': 8, 'wcs': None}
    args = auto_script.generate_fit_params(pline=pline, field_root=root, min_sens=0.0, min_mask=0.0)

    t0 = time.time()

    for src_id in sources:
        grizli_fit_z_niriss(root, src_id, fwhm, z1=z1, z2=z2)

    t1 = time.time()
    print("%.2f seconds" % (t1-t0))


def initializer(root):
    global grp
    grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), 
                            catalog='{0}-ir.cat.fits'.format(root), 
                            cpu_count=-1, sci_extn=1, pad=800)
    print("Number of FLTs:")
    print(len(grp.FLTs))

def grizli_fit_all_z_async_niriss(HOME_PATH='/Volumes/data3/grizli/',
                       root='my_data',
                       sources = [],
                       N_CPU='',
                       fwhm=325, mag_limit=30,z1=0.05,z2=12.0):

    #os.chdir(os.path.join(HOME_PATH + root, 'Extractions'))



    #global grp
    #grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), 
    #                        catalog='{0}-ir.cat.fits'.format(root), 
    #                        cpu_count=-1, sci_extn=1, pad=800)
    #print("Number of FLTs:")
    #print(len(grp.FLTs))


    #if not sources: 
    #    all_ids = [id for flt in grp.FLTs for id in flt.object_dispersers]
    #    sources = list(set(all_ids))
    #    sources.sort()

    print("Computing redshifts for %i sources" % (len(sources)))


    pline={'kernel': 'square', 'pixfrac': 0.5, 'pixscale': 0.04, 'size': 8, 'wcs': None}
    args = auto_script.generate_fit_params(pline=pline, field_root=root, min_sens=0.0, min_mask=0.0)

    # version 2 
    t0 = time.time()
    # Step 1: Init multiprocessing.Pool()
    print()
    print("Number of CPUs =",mp.cpu_count())
    if not N_CPU:
        N_CPU = mp.cpu_count()
    print("Using %i CPUs" % N_CPU)
    #pool = mp.Pool(mp.cpu_count())
    #pool = mp.Pool(N_CPU)
    pool = mp.Pool(N_CPU, initializer, (root,))

    processes = []
    # Step 3: Use loop to parallelize
    for src_id in sources:
        p = pool.apply_async(grizli_fit_z_niriss, args=(root,src_id,fwhm))
        processes.append(p)

    pool.close()
    pool.join()

    result = [p.get() for p in processes]

    t1 = time.time()
    print("%.2f seconds" % (t1-t0))

def grizli_write_all_beams(root='my_data',
                           sources = []):

    print("Writing beams for %i sources" % (len(sources)))


    grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), 
                            catalog='{0}-ir.cat.fits'.format(root), 
                            cpu_count=-1, sci_extn=1, pad=800)

    print("Number of FLTs:")
    print(len(grp.FLTs))

    for src_id in sources:
        print(src_id)

        try:

            beams = grp.get_beams(src_id, size=32, min_mask=0, min_sens=0)
            mb = multifit.MultiBeam(beams, fcontam=0.1, min_sens=0, min_mask=0, group_name=root)
            mb.write_master_fits()

        except:

            print("%s missing from model. Skipping." % (src_id))

            continue


def grizli_fit_z_niriss(root, id, fwhm, z1=0.05, z2=12.0):
    print(id)

    #beams = grp.get_beams(id, size=32, min_mask=0, min_sens=0)
    #mb = multifit.MultiBeam(beams, fcontam=0.1, min_sens=0, min_mask=0, group_name=root)
    #mb.write_master_fits()

    _fit = fitting.run_all_parallel(id, zr=[z1, z2], verbose=True, get_output_data=True)

    collected = gc.collect()
    print("Garbage collector: collected %d objects." % (collected))

    return 1
