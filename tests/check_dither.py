


import glob, sys
import astropy.io.fits as pyfits
import numpy as np


ext = 1

search_str = sys.argv[1]

L = glob.glob(search_str)
L.sort()
print(L)

ra = []
dec = []
for l in L:
    hdu = pyfits.open(l)
    #hdu.info()
    head = hdu[ext].header
    #print(head)
    #print(l)
    ra += [head["CRVAL1S"]]
    dec += [head["CRVAL2S"]]


# RA        Dec      delta_RA delta_Dec
#                    [arcsec]  [arcsec]                        


# old test
#ra = np.array([214.47307, 214.47507, 214.47307, 214.47107])
#dec = np.array([59.44919, 59.45119, 59.44769, 59.44619])
ind = np.arange(len(ra)) + 1

# new test
ra = np.array(ra)
dec = np.array(dec)
#dec -= 59.0
#dec = np.zeros((4,))
print(ra)
print(dec)
print()

d_ra = (ra - ra[0])*3600.*np.cos(dec[0]*np.pi/180.)
d_dec = (dec - dec[0])*3600.
print(d_ra)
print(d_dec)

N = len(ra)
print('%.8f %.8f %.2f" %.2f"' % (ra[0], dec[0], 0.0, 0.0))
for i in np.arange(1,N):
    d_ra = (ra[i] - ra[i-1])*3600.*np.cos(dec[i]*np.pi/180.)
    d_dec = (dec[i] - dec[i-1])*3600.

    print('%.8f %.8f %.2f" %.2f"' % (ra[i], dec[i], d_ra, d_dec))



radec_list = [ind, ra, dec]


def ds9_tmp(table_or_list,ra_key,dec_key,txt_key,ds9_file="catalogs/tmp.reg",color="green",rad='3"',txt_str="%.3f"):

    f = open(ds9_file,"w")
    #print "fk5"
    f.write("fk5\n")

    ra = table_or_list[ra_key]
    dec = table_or_list[dec_key]
    txt = table_or_list[txt_key]

    for a,d,m in zip(ra,dec,txt):
        #print 'circle(%.10f, %.10f, 1") # text={%.2f"}' % (a,d,m)
        #f.write('circle(%.10f, %.10f, %s) # text={%r} color=%s\n' % (a,d,rad,m,color))
        obj_str = 'circle(%.10f, %.10f, %s) # text={%s} color=%s\n' % (a,d,rad,txt_str % m,color)
        f.write(obj_str)
    f.close()



ds9_tmp(radec_list, 1, 2, 0, ds9_file="tmp.reg", rad='1"', txt_str="%i")
