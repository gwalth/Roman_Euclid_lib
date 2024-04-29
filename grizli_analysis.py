"""
Grizli Analysis
"""

import gc
import glob,os,time
import pickle

from collections import Counter, OrderedDict

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker

import numpy as np

import astropy.io.fits as pyfits
from astropy.table import Table, join, vstack
from astropy.coordinates import SkyCoord
import astropy.units as u

from grizli import multifit


emlines = [["OVI",          1038.0],         # 0
           ["Ly$\\alpha$",  1215.67],        # 1
           ["CIV",          1550.0],             # 2
           ["CIII]",        1909.],              # 3
           ["CII]",         2327.],              # 4
           ["MgII",         2796.4],             # 5
           ["MgII",         2803.5],             # 6
           ["NeV",          3326.],              # 7
           ["[OII]",        3727.],  # O2        # 8
           ["[NeIII]",      3868.7],             # 9
           ["H$\gamma$",    4340.5],  # Hg    # 10
           ["[OIII]",       4363.0],  # O31      # 11
           ["H$\\beta$",    4861.3],  # Hb   # 12
           ["[OIII]",       4959.0],  # O32      # 13
           ["[OIII]",       5007.0],  # O33      # 14
           ["[NII]",        6548.1],             # 15
           ["H$\\alpha$",   6562.8],  # Ha   # 16
           ["[NII]",        6583.0],             # 17
           ["[SII]",        6717.0],             # 18
           ["[SII]",        6731.0],             # 19
           ["[SIII]",       9069.0],             
           ["[SIII]",       9545.0],
           ["P$\\delta$",  10049.8],  # Pd   # 20
           ["P$\\gamma$",  10938.0],  # Pg   # 21
           ["P$\\beta$",   12818.1],  # Pb   # 22
           ["P$\\alpha$",  18750.1],  # Pa   # 23 
           ["Br$\\delta$", 19440.0],  # Br-d (wikipedia, not exact)
           ["Br$\\gamma$", 21660.0],  # Br-g (wikipedia, not exact)
           ["Br$\\beta$",  26250.0],  # Br-b (wikipedia, not exact)
           ["Br$\\alpha$", 40510.0],  # Br-a (wikipedia, not exact) 
          ]

# http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/hydspec.html
# http://articles.adsabs.harvard.edu//full/1934ApJ....80...19M/0000022.000.html


def display_grizli(root, id, w0=0.8, w1=1.7, labels=1, y0=-1, y1=-1, z_in=0, path="", lw=2, 
                   fontsize=8):
    
    f_full = '{0}_{1:05d}.full.fits'.format(root, id)
    #print(f_full)
    #full_hdu = pyfits.open(f_full)
    full_hdu = pyfits.open(os.path.join(path,f_full))
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
    #oned_hdu = pyfits.open(f_1d)
    oned_hdu = pyfits.open(os.path.join(path,f_1d))
    #print(oned_hdu[1].header)
    print(oned_hdu.info())
    grism = Table(oned_hdu['GR'].data)
    grism.show_in_notebook()
    print(grism.colnames)
    #print()
    
    
    
    f_2d = '{0}_{1:05d}.stack.fits'.format(root, id)
    #print(f_2d)
    #twod_hdu = pyfits.open(f_2d)
    twod_hdu = pyfits.open(os.path.join(path,f_2d))
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
    
    
    
    if z_in: z0 = z_in
    else: z0 = zfit_head['z_map']
        
    #z0 = zfit_head['z_map']
       
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
    multifit.show_drizzle_HDU(twod_hdu,scale_size=2.0)

def redshift_analysis_od(root,id,save_lines=3):

    full_obj = OrderedDict()
    
    fits_file = '{0}_{1:05d}.full.fits'.format(root, id)
   
    full_hdu = pyfits.open(fits_file)
    h0 = full_hdu[0].header
    line_keys = ['ID','RA','DEC','REDSHIFT','NUMLINES']
    #print(h0)
    for lk in line_keys:
        #h0[lk],h0.comments[lk]
        #print(h0[lk],)
        full_obj[lk] = h0[lk]
    
    # read all of the lines
    #if h0['NUMLINES'] < 1:
    #    print(root)
    #    print(h0['NUMLINES'])
    #    print(fits_file)
    #    return
        
    lines_list = []
    if h0['NUMLINES'] > 0:
        for i in range(h0['NUMLINES']):
            #print(i+1)
            od = OrderedDict({"line":h0["LINE%03d" % (i+1)],
                              "flux":h0["FLUX%03d" % (i+1)],
                              "err":h0["ERR%03d" % (i+1)]})                  
            lines_list.append(od)
        
    else:
        print(fits_file)
        print("NUMLINES =",h0['NUMLINES'])
        od = OrderedDict({"line": "None",
                          "flux": 99.9,
                          "err": -99.9})                  
        lines_list.append(od)
        
                                 
    #print(lines_list)
    line_tbl = Table(lines_list)
    line_tbl.sort("flux")
    line_tbl.reverse()
    
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


def plot_redshifts2(tbl,zcut=0.01, mag_key="mag_auto"):

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
    p2.scatter(tbl["REDSHIFT"], tbl[mag_key], s=2, c="k")
    p2.scatter(tbl_filt1["REDSHIFT"], tbl_filt1[mag_key])
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
    p4.scatter(tbl["REDSHIFT"], tbl[mag_key], s=2, c="k")
    p4.scatter(tbl_filt2["REDSHIFT"], tbl_filt2[mag_key])
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

##########
# checks #
##########

def check_primers(HOME_PATH = "/home/gwalth/data/Roman/grizli/sims",
                  LOCAL_PATH = "/local/RomanSims/grizli/sims",
                  root = "hlss", PREP_PATH = "Prep", glob_dir="", verb=0):
    
    
    print('HOME_PATH = ', HOME_PATH)
    print('LOCAL_PATH = ', LOCAL_PATH)
    
    comp = os.environ["HOST"]
    print(comp)

    #extract_path = os.path.join(LOCAL_PATH, root, 'Extractions/')
        
    os.chdir(os.path.join(HOME_PATH, root, PREP_PATH))


    if glob_dir: glob_dir += "/"
       
    # Primer catalog created for the aXeSIM 
    L = glob.glob(glob_dir+"Euclid_Roman_4deg2_*.fits")
    L.sort()
    primer_cats = [l for l in L if "slitless" not in l and "direct" not in l]
    print(primer_cats)
    
    for primer_cat in primer_cats:
        print(primer_cat)
        #primer = Table.read(primer_cat, format='ascii.sextractor')
        primer = Table.read(primer_cat, format='fits')
        #primer[:10].show_in_notebook()
        print(len(primer))
        print(primer.colnames)
        print()

    print(len(primer_cats))

def check_phot(HOME_PATH = "/home/gwalth/data/Roman/grizli/sims",
                   LOCAL_PATH = "/local/RomanSims/grizli/sims", mag_key = "mag_auto",
                   EXTRACT_PATH = "Extraction/", PREP_PATH = "Prep", cat = "",
                   prefix = "Euclid_Roman_4deg2_", suffix = "_final.cat",
                   root = "sim_v3", glob_dir="field*/", verb=0):
    
    
    print('HOME_PATH = ', HOME_PATH)
    print('LOCAL_PATH = ', LOCAL_PATH)
    
    comp = os.environ["HOST"]
    print(comp)

    extract_path = os.path.join(LOCAL_PATH, root, EXTRACT_PATH)
    extract_paths = glob.glob(extract_path + glob_dir)
    extract_paths.sort()
    #print(extract_paths)
    #print(len(extract_paths))
    #print()
    
    for extract_path in extract_paths:

        field = os.path.basename(extract_path)
        print("Processing %s" % field)
        
        os.chdir(os.path.join(HOME_PATH, root, PREP_PATH, field))
        
        if cat:
            ncat = cat
            format = "fits"
        else:
            # Photometry catalog created by SExtractor
            ncat = glob.glob(prefix + field + "*" + suffix)[0]
            format = 'ascii.sextractor'
        
        print(cat)
        phot = Table.read(ncat, format=format) # ref_cat in multimission
        print(len(phot))
        print(phot.colnames)
        print(len(phot.colnames))
        print()
        
        #print(np.min(phot["MAG_AUTO"]),np.max(phot["MAG_AUTO"]))

        if verb:
            fig = plt.figure()
            p = fig.add_subplot(111)
            p.hist(phot[mag_key],range=(10,32),bins=44)
            p.set_xlabel("H158 Magnitude")
            p.set_ylabel("N")
            p.set_xlim(10,32)

    print(len(extract_paths))

# checks for a single field
def check_single_phot(HOME_PATH = "/home/gwalth/data/Roman/grizli/sims",
                   LOCAL_PATH = "/local/RomanSims/grizli/sims",
                   cat = "hlss_phot.fits", mag_key = "mag_auto",
                   root = "hlss", verb=0):
    
    
    print('HOME_PATH = ', HOME_PATH)
    print('LOCAL_PATH = ', LOCAL_PATH)
    
    comp = os.environ["HOST"]
    print(comp)    
    os.chdir(os.path.join(HOME_PATH, root, 'Prep'))
       

    #phot = Table.read(cat, format='ascii.sextractor') # ref_cat in multimission
    phot = Table.read(cat, format='fits') 
    print(len(phot))
    print(phot.colnames)
    print(len(phot.colnames))
    print()
        
    #print(np.min(phot["MAG_AUTO"]),np.max(phot["MAG_AUTO"]))

    if verb:
        fig = plt.figure()
        p = fig.add_subplot(111)
        p.hist(phot[mag_key],range=(10,32),bins=44)
        p.set_xlabel("H158 Magnitude")
        p.set_ylabel("N")
        p.set_xlim(10,32)


def check_redshift_fits(HOME_PATH = "/home/gwalth/data/Roman/grizli/sims",
                        LOCAL_PATH = "/local/RomanSims/grizli/sims",
                        EXTRACT_PATH = "Extraction/", prefix = "sim_v3",
                        root = "sim_v3", glob_dir="field*", verb=0):

    
    TT0 = time.time()
    
    print('HOME_PATH = ', HOME_PATH)
    print('LOCAL_PATH = ', LOCAL_PATH)
    
    comp = os.environ["HOST"]
    print(comp)

    extract_path = os.path.join(LOCAL_PATH, root, EXTRACT_PATH)
    extract_paths = glob.glob(extract_path + glob_dir)
    extract_paths.sort()
    #print(extract_paths)
    #print(len(extract_paths))
    print()
    
    #for extract_path in extract_paths[6:10]: # problem extraction files in field_0017
    for extract_path in extract_paths:

        t0 = time.time()
        field = os.path.basename(extract_path)
        print("Processing %s" % field)
        
        os.chdir(os.path.join(LOCAL_PATH, root, EXTRACT_PATH, field))
        
        print("Find ids of all the extracted objects")
        ids = []
        files = glob.glob('*full.fits')
        for file in files:
            #print(file)
            ids.append(int(file.replace(".full.fits","").split("_")[-1]))

        #print(files)

        ids.sort()
        #print(ids)
        print("N =",len(ids))
        #print()

        print("Load redshift FITS files")
        redshift_rows = []
        for id in ids:
            redshift_od = redshift_analysis_od(prefix,id,save_lines=5) 
            redshift_od["FIELD"] = field
            redshift_rows.append(redshift_od)

        redshift_fits = Table(redshift_rows)
        print(redshift_fits.colnames)
        t1 = time.time()
        print("Time to load FITS files and create Table %.1f seconds" % (t1-t0))
        print()

    TT1 = time.time()
    
    print("Total time to process ALL files %.1f seconds" % (TT1 - TT0))
    print()
    print()

    print(len(extract_paths))

def check_roll_redshift_fits(HOME_PATH = "/home/gwalth/data/Roman/grizli/sims",
                             LOCAL_PATH = "/local/RomanSims/grizli/sims",
                             root = "hlss", prefix="", verb=0):
    
    if not prefix: prefix = root
    
    
    print('HOME_PATH = ', HOME_PATH)
    print('LOCAL_PATH = ', LOCAL_PATH)
    
    comp = os.environ["HOST"]
    print(comp)

    extract_path = os.path.join(LOCAL_PATH, root, 'Extractions/')
    os.chdir(os.path.join(LOCAL_PATH, root, 'Extractions'))

    t0 = time.time()

        
    print("Find ids of all the extracted objects")
    ids = []
    files = glob.glob('*full.fits')
    for file in files:
        #print(file)
        ids.append(int(file.replace(".full.fits","").split("_")[-1]))

    #print(files)

    ids.sort()
    #print(ids)
    print("N =",len(ids))
    #print()

    print("Load redshift FITS files")
    redshift_rows = []
    for id in ids:
        redshift_od = redshift_analysis_od(prefix,id,save_lines=5) 
        redshift_rows.append(redshift_od)

    redshift_fits = Table(redshift_rows)
    print(redshift_fits.colnames)
    t1 = time.time()
    print("Time to load FITS files and create Table %.1f seconds" % (t1-t0))
    print()

##########
# builds #
##########

def build_grp_tables(HOME_PATH = "/home/gwalth/data/Roman/grizli/sims",
                     LOCAL_PATH = "/local/RomanSims/grizli/sims",
                     EXTRACT_PATH = "Extractions/", PREP_PATH = "Prep", cat = "",
                     prefix = "Euclid_Roman_4deg2_", suffix = "_final.cat",
                     root = "sim_v3", mag_limit = 30, search_rad = 0.4,
                     glob_dir="field*", verb = 0):
    
    TT0 = time.time()
    
    print('HOME_PATH = ', HOME_PATH)
    print('LOCAL_PATH = ', LOCAL_PATH)
    
    comp = os.environ["HOST"]
    print(comp)

    extract_path = os.path.join(LOCAL_PATH, root, EXTRACT_PATH)
    extract_paths = glob.glob(extract_path + glob_dir)
    extract_paths.sort()
    print(extract_paths)
    #print(len(extract_paths))
    print()
    
    all_match_tables = []
    all_match_clean_tables = []
    
    #for extract_path in extract_paths[6:10]: # problem extraction files in field_0017
    for extract_path in extract_paths:
        T0 = time.time()
        
        field = os.path.basename(extract_path)
        print("Processing %s" % field)
        
        os.chdir(os.path.join(HOME_PATH, root, PREP_PATH, field))

        # Primer catalog created for the aXeSIM
        primer_cat = "Euclid_Roman_4deg2_%s.fits" % (field)
        print(primer_cat)
        #primer = Table.read(primer_cat, format='ascii.sextractor')
        primer = Table.read(primer_cat, format='fits')
        #primer[:10].show_in_notebook()
        print(len(primer))
        #print(primer.colnames) 
    
        if cat:
            ncat = cat
            format = "fits"
        else:
            # Photometry catalog created by SExtractor
            ncat = glob.glob(prefix + field + "*" + suffix)[0]
            format = 'ascii.sextractor'
        
        print(cat)
        phot = Table.read(ncat, format=format) # ref_cat in multimission
        print(len(phot))
        print()

        #print(np.min(phot["MAG_AUTO"]),np.max(phot["MAG_AUTO"]))

        if verb:
            fig = plt.figure()
            p = fig.add_subplot(111)
            p.hist(phot["MAG_AUTO"],range=(10,32),bins=44)
            p.set_xlabel("H158 Magnitude")
            p.set_ylabel("N")
            p.set_xlim(10,32)
              
        t0 = time.time()

        grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), 
                                catalog='{0}-ir.cat.fits'.format(prefix), 
                                cpu_count=-1, sci_extn=1, pad=800) 
    
        t1 = time.time()

        # all the ids in the grp (GLW)
        all_ids = [id for flt in grp.FLTs for id in flt.object_dispersers]
        print("Number of FLTs:")
        print(len(grp.FLTs))
        
        ids = list(set(all_ids))
        ids.sort()
        print("ids =")
        print(ids)
        print("Number of ids:")
        print(len(ids))
        
        sources = ids
        print(sources)
        print(len(sources))
        
        os.chdir(os.path.join(LOCAL_PATH, root, EXTRACT_PATH, field))
        bids = beams_ids(prefix)
        #print(bids)
        print(len(bids))
        help(beams_ids)
        
        #bids.index(595)
        sources = bids
        print(sources)
        print(len(sources))
        
        ############################################ 
            
        print("Find ids of all the extracted objects")
        ids = []
        files = glob.glob('*full.fits')
        for file in files:
            #print(file)
            ids.append(int(file.replace(".full.fits","").split("_")[-1]))

        #print(files)

        ids.sort()
        #print(ids)
        print("N =",len(ids))
        #print()

        
        t0 = time.time()
        print("Load redshift FITS files")
        redshift_rows = []
        for id in ids:
            redshift_od = redshift_analysis_od(prefix,id,save_lines=5) 
            redshift_od["FIELD"] = field
            redshift_rows.append(redshift_od)

        redshift_fits = Table(redshift_rows)
        print(redshift_fits[:10])
        print(phot[:10])
        print(len(redshift_fits))
        print(len(phot))

        t1 = time.time()

        print("Time to load FITS files and create Table %.1f seconds" % (t1-t0))
        print()
        
        print("Join primer and redshift_fits tables -> all_tbl")
        redshift_fits.rename_column('ID', 'id')
        all_tbl = join(phot, redshift_fits, keys='id')
        print(all_tbl)

        print("Search primer and all_tbl for RA/DEC matches")
        c_prime = SkyCoord(ra=primer["RA"]*u.degree, dec=primer["DEC"]*u.degree)
        c_all = SkyCoord(ra=all_tbl["RA"]*u.degree, dec=all_tbl["DEC"]*u.degree)

        idx, d2d, d3d = c_prime.match_to_catalog_sky(c_all)


        #print(idx)
        #print(d2d)
        #print(len(idx))
        #print(len(c_prime))
        #print(len(c_all))

        filt = d2d < 1*u.arcsec

        if verb:
            fig = plt.figure()
            p1 = fig.add_subplot(111)
            p1.hist(d2d[filt].value*3600.,bins=20)
        
        primer['idx'] = idx
        primer['d2d'] = d2d
        #print(primer.colnames)
        
        all_tbl['idx'] = np.arange(len(all_tbl))
        #print(all_tbl.colnames)
        #print(all_tbl['idx'])
        
        print("Join all_tbl and primer tables")
        match_tbl = join(all_tbl, primer, keys='idx')
        #print(match_tbl)
        #print(match_tbl.colnames)
        
        print("Select only sources < %.2f arcsec" % (search_rad))
        #filt = match_tbl['d2d'] < 1.0*u.arcsec
        clean_filt = match_tbl['d2d'] < search_rad*u.arcsec
        match_clean_tbl = match_tbl[clean_filt]
        print(len(match_clean_tbl))
        
        all_match_tables.append(match_tbl)
        all_match_clean_tables.append(match_clean_tbl)
        
        
        T1 = time.time()
        print()
        print("Total time to process %s = %.1f seconds" % (field,T1-T0))
        print()
        print()
        
    TT1 = time.time()
    
    print("Total time to process ALL files %.1f seconds" % (TT1 - TT0))
    print()
    print()
        
    return all_match_tables,all_match_clean_tables


def build_all_tables(HOME_PATH = "/home/gwalth/data/Roman/grizli/sims",
                     LOCAL_PATH = "/local/RomanSims/grizli/sims",
                     root = "sim_v3", mag_limit = 30, search_rad = 0.4,
                     glob_dir="field*", verb = 0):
    
    TT0 = time.time()
    
    print('HOME_PATH = ', HOME_PATH)
    print('LOCAL_PATH = ', LOCAL_PATH)
    
    comp = os.environ["HOST"]
    print(comp)

    extract_path = os.path.join(LOCAL_PATH, root, 'Extraction/')
    extract_paths = glob.glob(extract_path + glob_dir)
    extract_paths.sort()
    #print(extract_paths)
    #print(len(extract_paths))
    print()
    
    all_match_tables = []
    all_match_clean_tables = []
    
    #for extract_path in extract_paths[6:10]: # problem extraction files in field_0017
    for extract_path in extract_paths:
        T0 = time.time()
        
        field = os.path.basename(extract_path)
        print("Processing %s" % field)
        
        os.chdir(os.path.join(HOME_PATH, root, 'Prep', field))

        # Primer catalog created for the aXeSIM
        primer_cat = "Euclid_Roman_4deg2_%s.fits" % (field)
        print(primer_cat)
        #primer = Table.read(primer_cat, format='ascii.sextractor')
        primer = Table.read(primer_cat, format='fits')
        #primer[:10].show_in_notebook()
        print(len(primer))
        #print(primer.colnames) 
    
        prefix = new_direct.replace(".fits","")
        cat = prefix + ".cat" 
        print(cat)
        phot = Table.read(cat, format='ascii.sextractor') # ref_cat in multimission
        print(len(phot))
        
        
        #print(np.min(phot["MAG_AUTO"]),np.max(phot["MAG_AUTO"]))

        if verb:
            fig = plt.figure()
            p = fig.add_subplot(111)
            p.hist(phot["MAG_AUTO"],range=(10,32),bins=44)
            p.set_xlabel("H158 Magnitude")
            p.set_ylabel("N")
            p.set_xlim(10,32)
              
        t0 = time.time()

        print("Load Roman Grizli pickle")
        with open('Roman_GrismFLT.pickle', 'rb') as f:
            # The protocol version used is detected automatically, so we do not
            # have to specify it.
            Roman = pickle.load(f)
    
        t1 = time.time()

        print("Time to load pickle %.1f seconds" % (t1-t0))
        
        #print(Roman)
        #print(Roman.__dict__)
        #print(Roman.keys())

        #print()
        print("Check simulation")
        Roman_all,Roman_magcut,Roman_extract = check_sims(Roman, mag_limit)
        
        del Roman

        collected = gc.collect()
        print("Garbage collector: collected %d objects." % (collected))
        #print()
        
        os.chdir(os.path.join(LOCAL_PATH, root, 'Extraction', field))

        print("Find ids of all the extracted objects")
        ids = []
        files = glob.glob('*full.fits')
        for file in files:
            #print(file)
            ids.append(int(file.replace(".full.fits","").split("_")[-1]))

        #print(files)

        ids.sort()
        #print(ids)
        print("N =",len(ids))
        #print()

        
        t0 = time.time()
        print("Load redshift FITS files")
        redshift_rows = []
        for id in ids:
            redshift_od = redshift_analysis_od(root,id,save_lines=5)
            redshift_od["FIELD"] = field
            redshift_rows.append(redshift_od)

        redshift_fits = Table(redshift_rows)
        #print(redshift_fits)

        t1 = time.time()

        print("Time to load FITS files and create Table %.1f seconds" % (t1-t0))
        print()
        
        print("Join primer and redshift_fits tables -> all_tbl")
        Roman_extract.rename_column('NUMBER', 'ID')
        all_tbl = join(Roman_extract, redshift_fits, keys='ID')
        #print(all_tbl)

        print("Search primer and all_tbl for RA/DEC matches")
        c_prime = SkyCoord(ra=primer["RA"]*u.degree, dec=primer["DEC"]*u.degree)
        c_all = SkyCoord(ra=all_tbl["RA"]*u.degree, dec=all_tbl["DEC"]*u.degree)

        idx, d2d, d3d = c_prime.match_to_catalog_sky(c_all)


        #print(idx)
        #print(d2d)
        #print(len(idx))
        #print(len(c_prime))
        #print(len(c_all))

        filt = d2d < 1*u.arcsec

        if verb:
            fig = plt.figure()
            p1 = fig.add_subplot(111)
            p1.hist(d2d[filt].value*3600.,bins=20)
        
        primer['idx'] = idx
        primer['d2d'] = d2d
        #print(primer.colnames)
        
        all_tbl['idx'] = np.arange(len(all_tbl))
        #print(all_tbl.colnames)
        #print(all_tbl['idx'])
        
        print("Join all_tbl and primer tables")
        match_tbl = join(all_tbl, primer, keys='idx')
        #print(match_tbl)
        #print(match_tbl.colnames)
        
        print("Select only sources < %.2f arcsec" % (search_rad))
        #filt = match_tbl['d2d'] < 1.0*u.arcsec
        clean_filt = match_tbl['d2d'] < search_rad*u.arcsec
        match_clean_tbl = match_tbl[clean_filt]
        print(len(match_clean_tbl))
        
        all_match_tables.append(match_tbl)
        all_match_clean_tables.append(match_clean_tbl)
        
        
        T1 = time.time()
        print()
        print("Total time to process %s = %.1f seconds" % (field,T1-T0))
        print()
        print()
        
    TT1 = time.time()
    
    print("Total time to process ALL files %.1f seconds" % (TT1 - TT0))
    print()
    print()
        
    return all_match_tables,all_match_clean_tables

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

def beams_ids(root):
    """Generates list of ids from beams written
    
    Parameters
    ----------
    root: str
        Beam format 
            '{root}_{id}.beams.fits'
    
    Returns
    -------
    ids: list
    """
    
    L = glob.glob("{0}*.beams.fits".format(root))
    L.sort()
    
    ids = []
    for l in L:
        id = int(l.replace(".beams.fits","").replace(root + "_",""))
        ids.append(id)
    return ids


def build_roll_tables(HOME_PATH = "/home/gwalth/data/Roman/grizli/sims",
                      LOCAL_PATH = "/local/RomanSims/grizli/sims",
                      root = "hlss", prefix="", mag_limit = 30, search_rad = 0.4,
                      cat = "hlss_phot.fits", mag_key = "mag_auto",
                      verb = 0):
    
    if not prefix: prefix = root
    
    TT0 = time.time()
    
    print('HOME_PATH = ', HOME_PATH)
    print('LOCAL_PATH = ', LOCAL_PATH)
    
    comp = os.environ["HOST"]
    print(comp)

    extract_path = os.path.join(LOCAL_PATH, root, 'Extractions/')
    
    all_match_tables = []
    all_match_clean_tables = []
    

    T0 = time.time()
        

    ############################################         
    os.chdir(os.path.join(HOME_PATH, root, 'Prep'))
    
    # Primer catalog created for the aXeSIM 
    L = glob.glob("Euclid_Roman_4deg2_*.fits")
    primer_cats = [l for l in L if "slitless" not in l and "direct" not in l]
    print(primer_cats)
    
    for primer_cat in primer_cats:
        print(primer_cat)
        #primer = Table.read(primer_cat, format='ascii.sextractor')
        primer = Table.read(primer_cat, format='fits')
        #primer[:10].show_in_notebook()
        print(len(primer))
        print(primer.colnames)
        print()
       

    #phot = Table.read(cat, format='ascii.sextractor') # ref_cat in multimission
    phot = Table.read(cat, format='fits') 
    print(len(phot))
    print(phot.colnames)
    print(len(phot.colnames))
    print()
        
    #print(np.min(phot["MAG_AUTO"]),np.max(phot["MAG_AUTO"]))

    if verb:
        fig = plt.figure()
        p = fig.add_subplot(111)
        p.hist(phot[mag_key],range=(10,32),bins=44)
        p.set_xlabel("H158 Magnitude")
        p.set_ylabel("N")
        p.set_xlim(10,32)
        

              
    t0 = time.time()
   
    grp = multifit.GroupFLT(grism_files=glob.glob('*GrismFLT.fits'), 
                            catalog='{0}-ir.cat.fits'.format(prefix), 
                            cpu_count=-1, sci_extn=1, pad=800) 
    
    t1 = time.time()
    print("Time to load FLTs %.1f seconds" % (t1-t0))
        
        
    # all the ids in the grp (GLW)
    all_ids = [id for flt in grp.FLTs for id in flt.object_dispersers]
    print("Number of FLTs:")
    print(len(grp.FLTs))

    ids = list(set(all_ids))
    ids.sort()
    print("ids =")
    print(ids)
    print("Number of ids:")
    print(len(ids))
    
    sources = ids
    print(sources)
    print(len(sources))
    
    os.chdir(os.path.join(LOCAL_PATH, root, 'Extractions'))
    bids = beams_ids(prefix)
    #print(bids)
    print(len(bids))
    help(beams_ids)

    #bids.index(595)
    sources = bids
    print(sources)
    print(len(sources))
    
    ############################################ 
        
    #print(grp)
    #print(grp.__dict__)
    #print(grp.keys())

    #print()
    #print("Check simulation")


    #collected = gc.collect()
    #print("Garbage collector: collected %d objects." % (collected))
    ##print()
        
    os.chdir(os.path.join(LOCAL_PATH, root, 'Extractions'))

    print("Find ids of all the extracted objects")
    ids = []
    files = glob.glob('*full.fits')
    for file in files:
        #print(file)
        ids.append(int(file.replace(".full.fits","").split("_")[-1]))

    #print(files)

    ids.sort()
    #print(ids)
    print("N =",len(ids))
    #print()

        
    t0 = time.time()
    print("Load redshift FITS files")
    redshift_rows = []
    for id in ids:
        redshift_od = redshift_analysis_od(prefix,id,save_lines=5)
        redshift_rows.append(redshift_od)

    redshift_fits = Table(redshift_rows)
    print(redshift_fits[:10])
    print(phot[:10])
    print(len(redshift_fits))
    print(len(phot))

    t1 = time.time()

    print("Time to load FITS files and create Table %.1f seconds" % (t1-t0))
    print()
    

        
    print("Join primer and redshift_fits tables -> all_tbl")
    redshift_fits.rename_column('ID', 'id')
    all_tbl = join(phot, redshift_fits, keys='id')
    print(all_tbl)

    #sys.exit()
        
    #####################################
    #####################################
    #####################################
    
    
    
    print("Search primer and all_tbl for RA/DEC matches")
    c_prime = SkyCoord(ra=primer["RA"]*u.degree, dec=primer["DEC"]*u.degree)
    c_all = SkyCoord(ra=all_tbl["RA"]*u.degree, dec=all_tbl["DEC"]*u.degree)

    idx, d2d, d3d = c_prime.match_to_catalog_sky(c_all)


    #print(idx)
    #print(d2d)
    #print(len(idx))
    #print(len(c_prime))
    #print(len(c_all))

    filt = d2d < 1*u.arcsec

    if verb:
        fig = plt.figure()
        p1 = fig.add_subplot(111)
        p1.hist(d2d[filt].value*3600.,bins=20)
        
    primer['idx'] = idx
    primer['d2d'] = d2d
    #print(primer.colnames)
        
    all_tbl['idx'] = np.arange(len(all_tbl))
    #print(all_tbl.colnames)
    #print(all_tbl['idx'])
        
    print("Join all_tbl and primer tables")
    match_tbl = join(all_tbl, primer, keys='idx')
    #print(match_tbl)
    #print(match_tbl.colnames)
        
    print("Select only sources < %.2f arcsec" % (search_rad))
    #filt = match_tbl['d2d'] < 1.0*u.arcsec
    clean_filt = match_tbl['d2d'] < search_rad*u.arcsec
    match_clean_tbl = match_tbl[clean_filt]
    print(len(match_clean_tbl))
        
    #all_match_tables.append(match_tbl)
    #all_match_clean_tables.append(match_clean_tbl)
        
        
    TT1 = time.time()
    
    print("Total time to process ALL files %.1f seconds" % (TT1 - TT0))
    print()
    print()
        
    return match_tbl,match_clean_tbl

############
# analysis #
############

def sample_numbers(tbl, sn_key='SN001', sn=6.5, sigma_dz=0.005, zmin=0, zmax=6, alpha=1.0):
    
    filt = tbl[sn_key] > sn
    
    print('Num. sources =',len(tbl[sn_key]))
    print('Num. sources (S/N > %.1f) =' % (sn),len(tbl[sn_key][filt]))
    
    #z_fit = tbl["REDSHIFT_1"][filt]
    #z_true = tbl["REDSHIFT_2"][filt]
    z_fit = tbl["REDSHIFT"][filt]
    z_true = tbl["z_true"][filt]

    dz = (z_true - z_fit)/(1+z_true)

    #zfilt = np.abs(dz) < sigma_dz
    zfilt = np.abs(dz) < sigma_dz*(1+z_true)

    tbl_sn = tbl[filt]
    tbl_sn_dz = tbl_sn[zfilt]

    print("Num. sources (dz < %.3f) =" % (sigma_dz),len(z_fit[zfilt]))
    print("Frac. sources (dz < %.3f) = %.4f" % (sigma_dz,len(z_fit[zfilt])/len(z_fit)))
    print()
    
    fig = plt.figure(figsize=(4,4))
    p1 = fig.add_subplot(111)
    p1.hist(z_fit,bins=15,label="all",range=(zmin,zmax),alpha=alpha)
    p1.hist(z_fit[zfilt],bins=15,label="z < %.3f" % (sigma_dz),range=(zmin,zmax),alpha=alpha)
    p1.set_xlabel("Redshift")
    p1.set_ylabel("N")

    p1.legend()
    plt.show()
    
    return tbl_sn, tbl_sn_dz, filt, zfilt
    #return tbl_sn_dz
    

def sample_numbers2(tbl, sn_key='SN001', sn=6.5, sigma_dz=0.005, bin_size=0.1, zmin=0, zmax=6,
                    alpha=1.0, weights=None, full_output=0, label = "all", label_dz = "z < 0.005"):
    
    filt = tbl[sn_key] > sn
    
    print('Num. sources =',len(tbl[sn_key]))
    print('Num. sources (S/N > %.1f) =' % (sn),len(tbl[sn_key][filt]))
    
    #z_fit = tbl["REDSHIFT_1"][filt]
    #z_true = tbl["REDSHIFT_2"][filt]
    z_fit = tbl["REDSHIFT"][filt]
    z_true = tbl["z_true"][filt]

    dz = (z_true - z_fit)/(1+z_true)

    #zfilt = np.abs(dz) < sigma_dz
    zfilt = np.abs(dz) < sigma_dz*(1+z_true)

    tbl_sn = tbl[filt]
    tbl_sn_dz = tbl_sn[zfilt]

    print("Num. sources (dz < %.3f) =" % (sigma_dz),len(z_fit[zfilt]))
    print("Frac. sources (dz < %.3f) = %.4f" % (sigma_dz,len(z_fit[zfilt])/len(z_fit)))
    print()


    bins = int((zmax-zmin)/bin_size)
    print(bins)
    
    fig = plt.figure(figsize=(4,4))
    p1 = fig.add_subplot(111)


    hist, bin_edges = np.histogram(z_fit, bins=bins, range=(zmin,zmax))

    hist_dz, bin_edges_dz = np.histogram(z_fit[zfilt], bins=bins, range=(zmin,zmax))


    #p1.bar(bin_edges[:-1], hist, align="edge", width=bin_size)
    #p1.bar(bin_edges_dz[:-1], hist_dz, align="edge", width=bin_size)

    p1.bar(bin_edges_dz[:-1], hist_dz/hist, align="edge", width=bin_size)

    print(hist)
    print(hist_dz)
    print(bin_edges)
    print(bin_edges_dz)


    p1.set_xlabel("Redshift")
    p1.set_ylabel("N")

    p1.set_ylim(0,1.1)
    p1.set_xlim(zmin,zmax)

    p1.legend()
    plt.show()

   
    return tbl_sn, tbl_sn_dz, filt, zfilt
    
def sample_numbers3(tbl, sn_key='SN001', sn=6.5, sigma_dz=0.005, bin_size=0.1, zmin=0, zmax=6,
                    alpha=1.0, weights=None, full_output=0, label = "all", label_dz = "z < 0.005",
                    ax=""):
    
    filt = tbl[sn_key] > sn
    
    print('Num. sources =',len(tbl[sn_key]))
    print('Num. sources (S/N > %.1f) =' % (sn),len(tbl[sn_key][filt]))
    
    #z_fit = tbl["REDSHIFT_1"][filt]
    #z_true = tbl["REDSHIFT_2"][filt]
    z_fit = tbl["REDSHIFT"][filt]
    z_true = tbl["z_true"][filt]

    dz = (z_true - z_fit)/(1+z_true)

    #zfilt = np.abs(dz) < sigma_dz
    zfilt = np.abs(dz) < sigma_dz*(1+z_true)

    tbl_sn = tbl[filt]
    tbl_sn_dz = tbl_sn[zfilt]

    print("Num. sources (dz < %.3f) =" % (sigma_dz),len(z_fit[zfilt]))
    print("Frac. sources (dz < %.3f) = %.4f" % (sigma_dz,len(z_fit[zfilt])/len(z_fit)))
    print()


    bins = int((zmax-zmin)/bin_size)
    print(bins)

    plt_show = 0
    if not ax:
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111)
        plt_show = 1


    hist, bin_edges = np.histogram(z_fit, bins=bins, range=(zmin,zmax))

    hist_dz, bin_edges_dz = np.histogram(z_fit[zfilt], bins=bins, range=(zmin,zmax))

    bin_mid = [(bin_edges[i]+bin_edges[i+1])/2. for i in range(len(bin_edges[:-1]))]
    bin_mid_dz = [(bin_edges_dz[i]+bin_edges_dz[i+1])/2. for i in range(len(bin_edges_dz[:-1]))]


    #ax.bar(bin_edges[:-1], hist, align="edge", width=bin_size)
    #ax.bar(bin_edges_dz[:-1], hist_dz, align="edge", width=bin_size)

    #ax.bar(bin_edges_dz[:-1], hist_dz/hist, align="edge", width=bin_size)
    ax.plot(bin_mid_dz, hist_dz/hist, drawstyle="steps-mid",c="orange")
    #ax.plot(bin_edges_dz[:-1], hist_dz/hist, drawstyle="steps-pre",c="g")
    #ax.plot(bin_edges_dz[:-1], hist_dz/hist, drawstyle="steps-post",c="y")
    #ax.plot(bin_edges_dz[:-1], hist_dz/hist, drawstyle="steps-mid",c="r")

    print(hist)
    print(hist_dz)
    print(bin_edges)
    print(bin_edges_dz)

    ax.plot([zmin,zmax],[0.5,0.5],"--",lw=0.5,c="k")
    ax.plot([zmin,zmax],[0.9,0.9],"-.",lw=0.5,c="k")


    ax.set_xlabel("Redshift z")
    ax.set_ylabel("Frac")

    ax.set_ylim(0,1.1)
    ax.set_xlim(zmin,zmax)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(bin_size))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))

    #ax.legend()

    if plt_show:
        plt.show()

   
    return tbl_sn, tbl_sn_dz, filt, zfilt
    




def emline_bar_chart(data,title=""):

    c1 = Counter(data)
    c2 = OrderedDict(c1.most_common())

    ind = np.arange(len(c1))
    width=0.75

    fig = plt.figure()
    p1 = fig.add_subplot(111)
    p1.bar(ind, c2.values(), width)
    p1.set_xticks(ind)
    p1.set_xticklabels(c2.keys(), rotation = 90)
    p1.set_title(title)
    plt.show()
    
def emline_bar_chart2(data1,data2,title=""):

    c1 = Counter(data1)
    c1o = OrderedDict(c1.most_common())
    
    c2 = Counter(data2)
    c2o = OrderedDict(c2.most_common())
    
   
    ind1 = np.arange(len(c1))
    ind2 = np.arange(len(c2))
    width=0.75

    fig = plt.figure(figsize=(12,6))
    p1 = fig.add_subplot(121)
    p1.bar(ind1, c1o.values(), width)
    p1.set_xticks(ind1)
    p1.set_xticklabels(c1o.keys(), rotation = 90)
    p1.set_title(title)
    
    p2 = fig.add_subplot(122)
    p2.bar(ind2, c2o.values(), width)
    p2.set_xticks(ind2)
    p2.set_xticklabels(c2o.keys(), rotation = 90)
    p2.set_title(title)
    plt.show()
    
    plt.show()

def sn_dz_plot(tbl, sn_key='SN001'):    
    
    z_fit = tbl["REDSHIFT"]
    z_true = tbl["z_true"]

    dz = (z_true - z_fit)/(1+z_true)
    
    sn = tbl[sn_key]
    
    
    fig = plt.figure(figsize=(14,4))
    p1 = fig.add_subplot(111)
    p1.scatter(dz,sn,s=1,c="k")
    #p1.plot([-4,1],[5,5],"--",c="red")
    p1.set_xlabel("dz")
    p1.set_ylabel("S/N")
    
    #p1.set_yscale("log")
    
    p1.set_xlim(-0.01,0.01)
    p1.set_ylim(-1,100)
    #p1.legend()
    plt.show()


def sn_dz_hist2d_plot(tbl, sn_key='SN001'):    
    
    z_fit = tbl["REDSHIFT"]
    z_true = tbl["z_true"]

    dz = (z_true - z_fit)/(1+z_true)
    
    sn = tbl[sn_key]
    
    
    fig = plt.figure(figsize=(18,4))
    p1 = fig.add_subplot(111)
    im1 = p1.hist2d(dz,sn,bins=(75,20),range=[[-0.01,0.01],[0,30]],cmap="PuBu")
                    #cmap="cividis_r")
                    #cmap="magma_r")
                    #cmap="inferno_r")
                    #cmap="plasma_r")
                    #cmap="viridis_r")
    #p1.plot([-4,1],[5,5],"--",c="red")
    #p1.set_xlabel("dz")
    #p1.set_ylabel("S/N")
    
    #p1.set_yscale("log")
    
    #p1.set_xlim(-0.2,0.2)
    #p1.set_ylim(-1,10)
    #p1.legend()
    
    cb1 = fig.colorbar(im1[3])
    
    plt.show()
    
