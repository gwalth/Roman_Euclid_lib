#!/usr/bin/env python

import astropy.io.fits as pyfits
from astropy import wcs

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

filename = "/Users/gwalth/data/Roman/grizli/sims/Euclid/Env9_V1_TEST/Prep/NISP_SLITLESS_FRAME1_Env9_V1_TESTarray_noMW_RNDC_2024-02-05.fits"



# Load the FITS hdulist using astropy.io.fits
hdulist = pyfits.open(filename)
hdulist.info()


print()
print(hdulist[0].header)
print()
print(hdulist[1].header)
print()
print()
print()
print()

newhead1 = euclid_wcs(hdulist[1].header)

#print()
#print(newhead1)
#print()
#print(type(newhead1))
#print()
#print(dir(newhead1))
#print()
#print(newhead1.__dict__)
#print()


hdulist[0].header.update(newhead1)

remove_basic_wcs(hdulist[1].header, alt='S')
hdulist[1].header.update(newhead1)




print(hdulist[0].header)
print(hdulist[1].header)


#newhead = head0+newhead1
#newhead = head1+newhead1
#print(newhead)
