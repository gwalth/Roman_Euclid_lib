
import glob, os, time, sys
sys.path.append('/Users/gwalth/python/src')
import axe_disperse
from axe_disperse import aXeConf

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator

from astropy.table import Table, unique, join
import astropy.io.fits as pyfits
#from astropy import wcs
from astropy.wcs import WCS

from scipy import integrate, interpolate
from grizli_functions import euclid_wcs, remove_basic_wcs, stats, fit_gauss, fit_poly

c = 2.9979e10  # cm/s
Ang = 1e-8  # cm



def get_range(img, sigma=5.): 

    X = img.flatten()
    stat_dict = stats(X)
    vmin = stat_dict["median"] - sigma * stat_dict["mad"]
    vmax = stat_dict["median"] + sigma * stat_dict["mad"]

    return vmin, vmax


class aXeConf_glw(aXeConf):
    """extention of aXeConf for Roman and Euclid use"""


    #def extract1d(self, xc, yc, aper=10):
    
    def load_sensitivity(self, sens_file):


        self.sens_cat = Table.read(sens_file)
        #print(self.sens_cat.colnames)

        self.interp_sens = interpolate.interp1d(self.sens_cat['WAVELENGTH'],
                                                self.sens_cat['SENSITIVITY'])
        #w0, w1 = sens_cat['WAVELENGTH'][0], sens_cat['WAVELENGTH'][-1]

    def get_beam(self, xc, yc, beam='A', dx0 = -100, dx1 = 250):
        """
        Horizontal beam
        similar to aXe_conf.show_beams
        """
    
        dx = np.arange(dx0, dx1)

        xoff = self.field_dependent(xc, yc, self.conf['XOFF_%s' % (beam)])
        dy, lam = self.get_beam_trace(xc, yc, dx=dx, beam=beam)
        xlim = self.conf['BEAM%s' % (beam)]
        ok = (dx >= xlim[0]) & (dx <= xlim[1])

        #print(xoff)
        #print(dy)
        #print(lam)
        #print(xlim)
        #print(dy[ok])
        #print(lam[ok])

        x = xc + dx
        y = yc + dy
        #print(x)
        #print(y)
        
        self.x = x
        self.y = y
        self.xc = xc
        self.yc = yc
        self.dx = dx
        self.dy = dy
        self.lam = lam
        
        return x, y, dx, dy, lam
    
    def get_beam_invert(self, xc, yc, beam='A', dy0 = -550, dy1 = 550):
        """
        Vertical beam
        similar to aXe_conf.show_beams
        """
    
        dy = np.arange(dy0, dy1)

        yoff = self.field_dependent(xc, yc, self.conf['YOFF_%s' % (beam)])
        dx, lam = self.get_beam_trace(xc, yc, dx=dy, beam=beam)
        ylim = self.conf['BEAM%s' % (beam)]
        ok = (dy >= ylim[0]) & (dy <= ylim[1])

        #print(yoff)
        #print(dy)
        #print(lam)
        #print(ylim)
        #print(dx[ok])
        #print(lam[ok])

        x = xc + dx
        y = yc + dy
        #print(x)
        #print(y)
        
        self.x = x
        self.y = y
        self.xc = xc
        self.yc = yc
        self.dx = dx
        self.dy = dy
        self.lam = lam
        
        return x, y, dx, dy, lam
        
    def load_files(self, direct_file, slitless_file, pre_process_path=None, offset=0.0, fix_mosa=1, det='11',
remove_dc=0, correct_gain=0, gain=None, diag=1, exptime=574., clip_bad=1):

        self.offset = float(offset)
        self.exptime = float(exptime)

        if (pre_process_path is not None) and remove_dc:
            dc_file = os.path.join(pre_process_path,"DC_%s.fits" % det)
            print(dc_file)
            hdu = pyfits.open(dc_file)
            self.DC = hdu[0].data

        else:
            self.DC = 0.

        print("dc = ",self.DC)


        #if (pre_process_path is not None) and correct_gain:
        if (pre_process_path is not None) and correct_gain:

            if gain == None:
                gain_file = os.path.join(pre_process_path,"GAIN_%s.fits" % det)
                hdu = pyfits.open(gain_file)
                self.GAIN = np.mean(hdu[0].data)

            else:
                self.GAIN = gain

        else:
            self.GAIN = 1.0

        print("gain = ",self.GAIN)


        hdu = pyfits.open(slitless_file)
        


        raw = hdu['DET%s.SCI' % (det)].data
        chi2 = hdu['DET%s.CHI2' % (det)].data

        sci = (raw - self.offset)/self.exptime/self.GAIN
        #err = np.sqrt(chi2)/self.exptime/self.GAIN
        err = chi2/self.exptime/self.GAIN

        if remove_dc:
            sci -= self.DC


        if clip_bad:
            X = sci.flatten()
            Y = err.flatten()

            print("N =", len(X))
 
            stat_dict = stats(X)

            sigma = 10.
            x_med = stat_dict["median"]
            x_mean = stat_dict["mean"]
            x_clip_min = stat_dict["median"] - sigma * stat_dict["mad"]
            x_clip_max = stat_dict["median"] + sigma * stat_dict["mad"]
            print(x_clip_min, x_clip_max)

            filt = (X > x_clip_min) & (X < x_clip_max) 

            X_clip = X[filt]

            print("N =", len(X_clip))
            print("pixels clipped =", len(X) - len(X_clip))
            print("fraction of good pixels =", len(X_clip)/len(X))
             
            fig = plt.figure(figsize=(10,5))
          
            ax1 = fig.add_subplot(121)
            #ax1.hist(X, bins=20, range=(0,2000)) # Counts
            #ax1.hist(X, bins=200)
            ax1.hist(X, bins=100, range=(-1, 1))
            y0, y1 = ax1.get_ylim()

            ax1.plot([x_mean, x_mean], [y0,y1], c="tab:green", lw=2, label=r"$\bar{x}$") 
            ax1.plot([x_med, x_med], [y0, y1],"-", lw=2, c="tab:orange", label=r"$x_{\rm median}$")
            ax1.plot([x_clip_max, x_clip_max],[y0,y1],"--",lw=2, c="tab:orange",
                label=r"$x_{\rm median}+\sigma_{\rm MAD}$") 
            ax1.plot([x_clip_min, x_clip_min],[y0,y1],"--",lw=2, c="tab:orange",
                label=r"$x_{\rm median}-\sigma_{\rm MAD}$") 
            print(y0,y1)


            ax1.set_yscale("log")
            ax1.set_ylim(y0,y1)
            #ax1.set_xlabel("Counts")
            ax1.set_xlabel("ADU/s")
            ax1.set_title("SCI")
            #ax1.xaxis.set_major_locator(MultipleLocator(5))
          
            ax_ins = inset_axes(
                ax1, 
                width="100%", 
                height="100%",
                bbox_to_anchor=(0.125, 0.675, 0.3, 0.3), # right side
                bbox_transform=ax1.transAxes,
                borderpad=0
            )
          
            ax_ins.hist(X, bins=20)
            ax_ins.set_yscale("log")
            #ax_ins.set_xlabel("Counts")
            ax_ins.set_xlabel("ADU/s")
          
            ax_ins.xaxis.set_major_locator(MultipleLocator(10))
            #ax_ins.xaxis.set_minor_locator(MultipleLocator(0.1))
            #ax_ins.yaxis.set_major_locator(MultipleLocator(200))
            #ax_ins.yaxis.set_minor_locator(MultipleLocator(50))
          
            #ax_ins.tick_params(which='major', width=2)
            #ax_ins.tick_params(which='minor', width=1)

            leg = ax1.legend()
          
          
            ax2 = fig.add_subplot(122)
            #ax2.hist(Y,bins=20, range=(0,100)) # Counts
            ax2.hist(Y, bins=100)
            ax2.set_yscale("log")
            #ax2.set_xlabel("Counts")
            ax2.set_xlabel("ADU/s")
            ax2.set_title("ERR")
          
            plt.show()




            fig = plt.figure()
            
            print()
            vmin, vmax = get_range(sci, sigma=5.) 
            print(vmin, vmax)
            print()
            
            ax1 = fig.add_subplot(121)
            ax1.imshow(sci, vmin=vmin, vmax=vmax, origin="lower")
            
            #sci_clip = sci[sci < x_clip_min] = -99
            print("N(F < %.2f) =" % (x_clip_min), len(sci[sci < x_clip_min]))
            print("N(F < -0.25) =", len(sci[sci < -0.25]))
            print("N(F < -0.50) =", len(sci[sci < -0.50]))
            
            sci[sci < x_clip_min] = np.nan
            
            ax2 = fig.add_subplot(122)
            ax2.imshow(sci, vmin=vmin, vmax=vmax, origin="lower")
            
            plt.show()



        # sci, err
        self.slitless = [sci, err]
        self.head_slit = hdu['DET%s.SCI' % (det)].header


        if diag:

           fig = plt.figure()

           print()
           vmin, vmax = get_range(raw, sigma=5.) 
           print(vmin, vmax)
           print()

           ax1 = fig.add_subplot(121)
           ax1.imshow(raw, vmin=vmin, vmax=vmax, origin="lower")


           print()
           vmin, vmax = get_range(sci, sigma=5.) 
           print(vmin, vmax)
           print()

           ax2 = fig.add_subplot(122)
           ax2.imshow(sci, vmin=vmin, vmax=vmax, origin="lower")


           plt.show()

           #sys.exit()
           



        if fix_mosa:
            self.head_wcs = euclid_wcs( self.head_slit )
        else:
            self.head_wcs = self.head_slit 

        self.wcs_slit = WCS(self.head_wcs)
        print(self.wcs_slit)
        self.footprint_slit = self.wcs_slit.calc_footprint(header=self.head_slit)

        hdu = pyfits.open(direct_file)
        self.direct = [hdu[1].data, hdu[2].data, hdu[3].data]
        self.direct_slitless = hdu[1].header
        self.head_direct = hdu[1].header
        self.wcs_direct = WCS(header=self.head_direct)      
        self.footprint_direct = self.wcs_direct.calc_footprint(header=self.head_direct)

    def get_beam_cutout(self, xc, yc, aper, beam = 'A', dx0 = -100, dx1 = 250, offset=0.0):
        print("get_beam_cutout")
        
        x, y, dx, dy, lam = self.get_beam(xc, yc, beam=beam, dx0 = dx0, dx1 = dx1)
        
        self.aper = aper
        self.aper_hw = aper/2.

        print("y =", y)

        # the ends are causing problems with the interpolation 
        x0 = int(np.min(x) + 1) 
        x1 = int(np.max(x) - 1)
        y0 = int(np.min(y - self.aper_hw) + 0.5)
        y1 = int(np.max(y + self.aper_hw) + 0.5)
        
        interp_wav = interpolate.interp1d(x, lam)
        interp_y_dist = interpolate.interp1d(lam, y)

        print("(x0, x1) =", x0, x1)
        print("(y0, y1) =", y0, y1)
        print("(xc, yc) =", xc, yc)

        if x0 < 0: x0 = 0
        if x1 >= 2048: x1 = 2048

        # converting to prime coordinates (cutout coordinates)
        self.xp = xc - x0
        self.yp = yc - y0
        self.dxp = self.dx + self.xp
        self.dyp = self.dy + self.yp
        
        self.pix = np.arange(x0,x1,1)
        self.wav = [interp_wav(p) for p in self.pix]
        self.y_dist = np.array([interp_y_dist(w) for w in self.wav])



        self.cutout = [
            self.slitless[0][y0:y1,x0:x1] + offset, 
            self.slitless[1][y0:y1,x0:x1], 
        ]
        
        return self.cutout

    def get_beam_invert_cutout(self, xc, yc, aper, beam = 'A', dy0 = -100, dy1 = 250, offset=0.0):

        print("get_beam_invert_cutout")
        
        x, y, dx, dy, lam = self.get_beam_invert(xc, yc, beam=beam, dy0 = dy0, dy1 = dy1)
        
        self.aper = aper
        self.aper_hw = aper/2.
        #print("y =", y)

        # the ends are causing problems with the interpolation 
        y0 = int(np.min(y) + 1) 
        y1 = int(np.max(y) - 1)
        x0 = int(np.min(x - self.aper_hw) + 0.5)
        x1 = int(np.max(x + self.aper_hw) + 0.5)
        
        interp_wav = interpolate.interp1d(y, lam)
        interp_y_dist = interpolate.interp1d(lam, x)

        print("(x0, x1) =", x0, x1)
        print("(y0, y1) =", y0, y1)
        print("(xc, yc) =", xc, yc)

        if y0 < 0: y0 = 0
        if y1 >= 2048: y1 = 2048

        # converting to prime coordinates (cutout coordinates)
        self.xp = xc - x0
        self.yp = yc - y0
        self.dxp = self.dx + self.xp
        self.dyp = self.dy + self.yp
        
        self.pix = np.arange(y0,y1,1)
        self.wav = np.array([interp_wav(p) for p in self.pix])
        self.y_dist = np.array([interp_y_dist(w) for w in self.wav])


        self.cutout = [
            self.slitless[0][y0:y1,x0:x1] + offset, 
            self.slitless[1][y0:y1,x0:x1], 
        ]
        
        return self.cutout
    
    def get_thumbnail(self, a0, d0, xhw=50, yhw=50):
                
        p0 = [a0, d0]
        pixcoords = self.wcs_direct.wcs_world2pix([p0],1)
        xc = pixcoords[0][0]
        yc = pixcoords[0][1]
        
        yl,xl = self.direct[0].shape
        
        y0 = int(yc-yhw+0.5)
        y1 = int(yc+yhw+0.5)
        x0 = int(xc-xhw+0.5)
        x1 = int(xc+xhw+0.5)
        
        if y0 < 0: y0 = 0
        if y1 > yl: y1 = yl
        if x0 < 0: x0 = 0
        if x1 > xl: x1 = xl
        
        thumbnail = self.direct[0][y0:y1, x0:x1]
        extra = (x0, x1, y0, y1, xl, yl)
                
        self.thumbnail = thumbnail, extra

    
    def plot_thumbnail(self, cat=None, ra_key="X_WORLD", dec_key="Y_WORLD", 
                       label_key="MAG_AUTO", plot_labels=1, vmin=-1, vmax=1, verb=0):
        
        fig = plt.figure()    
        
        #if plot_all:
        #    print(labels.shape)
        #    filt = np.ones(labels.shape[0], dtype=int)
        #else:
        #    filt = labels < mag_cut
        
        thumbnail, extra = self.thumbnail
        x0, x1, y0, y1, xl, yl = extra

        ax1 = fig.add_subplot(111)
        ax1.imshow(thumbnail, origin="lower", vmin=vmin, vmax=vmax)
        
        if cat is not None:
            ra = cat[ra_key]
            dec = cat[dec_key]
            labels = cat[label_key]
            
            skycoords = np.array([[r0,d0] for r0,d0 in zip(ra,dec)])          
            pixcoords = self.wcs_direct.wcs_world2pix(skycoords, 0)

            #print("I am here")
            #print(pixcoords)
            
            filt = (pixcoords[:,0] > x0) & (pixcoords[:,0] < x1) & \
                   (pixcoords[:,1] > y0) & (pixcoords[:,1] < y1)
            
            #print(pixcoords[filt])
            #print(len(pixcoords[filt]))
            #print()
            
            x = pixcoords[:,0] - x0
            y = pixcoords[:,1] - y0
            
            #print(x[filt])
            #print(y[filt])
            #print(labels[filt])
        
            ax1.scatter(x[filt], y[filt], fc="None", ec="r", s=60, lw=0.3)
            
            if plot_labels:
                for j,label in enumerate(labels[filt]):
                    ax1.text(x[filt][j], y[filt][j], label)
                    
        plt.show()
    
    def plot_cutout(self, vmin=-1, vmax=1, verb=0, flip_axes=0):
        
        if verb:
            print("min =", np.nanmin(self.cutout[0]))
            print("max =", np.nanmax(self.cutout[0]))
            print("mean =", np.nanmean(self.cutout[0]))
            print("median =", np.nanmedian(self.cutout[0]))
            print("std =", np.nanstd(self.cutout[0]))
            print("q( 2.3) =", np.quantile(self.cutout[0], 0.023))
            print("q(15.9) =", np.quantile(self.cutout[0], 0.159))
            print("q(84.1) =", np.quantile(self.cutout[0], 0.841))
            print("q(97.7) =", np.quantile(self.cutout[0], 0.977))        
        
        #fig = plt.figure(figsize=(10,2.5))
        fig = plt.figure()

        print(self.x)
        print(self.y)
        print(self.dx)
        print(self.dy)
        print(self.xc)
        print(self.yc)

        ax1 = fig.add_subplot(211)
        if flip_axes:
            ax1.imshow(self.cutout[0].T, origin="lower", vmin=vmin, vmax=vmax)
            ax1.plot(self.dyp, self.dxp - 1, "--", c="tab:red")
            ax1.plot(self.dyp, self.dxp - self.aper_hw - 1, "-", c="tab:red")
            ax1.plot(self.dyp, self.dxp + self.aper_hw - 1, "-", c="tab:red")
            #ax1.plot(self.dy+self.yp, self.dx+self.xp-1, "--", c="tab:red")
            #ax1.plot(self.dy+self.yp, self.dx+self.xp-self.aper_hw-1., "-", c="tab:red")
            #ax1.plot(self.dy+self.yp, self.dx+self.xp+self.aper_hw-1., "-", c="tab:red")

            #y = x + c

            #filt1 = 
            
        else:
            ax1.imshow(self.cutout[0], origin="lower", vmin=vmin, vmax=vmax)
            ax1.plot(self.dxp, self.dyp - 1, "--", c="tab:red")
            ax1.plot(self.dxp, self.dyp - self.aper_hw - 1, "-", c="tab:red")
            ax1.plot(self.dxp, self.dyp + self.aper_hw - 1, "-", c="tab:red")
            #ax1.plot(self.dx+self.xp, self.dy+self.yp-1, "--", c="tab:red")
            #ax1.plot(self.dx+self.xp, self.dy+self.yp-self.aper_hw-1., "-", c="tab:red")
            #ax1.plot(self.dx+self.xp, self.dy+self.yp+self.aper_hw-1., "-", c="tab:red")

        ax1.set_aspect(3.)

        ax2 = fig.add_subplot(212)
        #ax2.plot(self.profile)
        ax2.plot(self.profile, drawstyle="steps-mid")
        #ax2.set_aspect("auto")

        plt.show()
        
    def plot_extraction(self, w0=None, w1=None, y0=None, y1=None, flux_units=False, model=None, scale=1.0, wav_units=True ):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        if flux_units:
            temp_flux = self.spectrum[2] * scale
            temp_err = self.spectrum[3] * scale
        else:
            temp_flux = self.spectrum[0] * scale
            temp_err = self.spectrum[1] * scale

        if wav_units:
            ax.plot(self.wav, temp_flux, label="Spectrum")
            ax.plot(self.wav, temp_err, label="Error")

        else:
            ax.plot(temp_flux, label="Spectrum")
            ax.plot(temp_err, label="Error")

        
        w0_lim, w1_lim = ax.get_xlim()
        y0_lim, y1_lim = ax.get_ylim()    
        
        if model != None and wav_units:
            wav, flux = model
            ax.plot(wav, flux, label="Model")
        
        if w0 == None: w0 = w0_lim
        if w1 == None: w1 = w1_lim
        if y0 == None: y0 = y0_lim
        if y1 == None: y1 = y1_lim

        if wav_units:            
            ax.set_xlim(w0,w1)
        ax.set_ylim(y0,y1)

        ax.legend()
        #ax1.plot(lam, spec)

        plt.show()
        
    def plot_trace(self, x0=None, x1=None, y0=None, y1=None, vmin=-0.05, vmax=0.05):

        fig = plt.figure()   
        ax1 = fig.add_subplot(111)
        ax1.imshow(self.slitless[0], origin="lower", vmin=vmin, vmax=vmax)

        #ax1.plot(x,y,"--",c="tab:red")
        print(self.xc,self.yc)
        ax1.plot(self.x, self.y - self.aper_hw, "--", c="tab:red")
        ax1.plot(self.x, self.y + self.aper_hw, "--", c="tab:red")
        
        x0_lim, x1_lim = ax1.get_xlim()
        y0_lim, y1_lim = ax1.get_ylim()       
        
        if x0 == None: x0 = x0_lim
        if x1 == None: x1 = x1_lim
        if y0 == None: y0 = y0_lim
        if y1 == None: y1 = y1_lim
            
        ax1.set_xlim(x0,x1)
        ax1.set_ylim(y0,y1)

        plt.show()
    
    def plot_footprint(self):
        self.footprint_slit = self.wcs_slit.calc_footprint(header=self.head_slit)
        self.footprint_direct = self.wcs_direct.calc_footprint(header=self.head_direct)

        X0,Y0 = 9999999,9999999
        X1,Y1 = -9999999,-9999999
        
        fig = plt.figure()   
        ax1 = fig.add_subplot(111)
        
        colors = ["tab:blue","tab:green"]
        
        for i, footprint in enumerate([self.footprint_direct, self.footprint_slit]):
        
            x0,y0 = np.min(footprint,axis=0)
            x1,y1 = np.max(footprint,axis=0)
    
            if x0 < X0: X0 = x0
            if y0 < Y0: Y0 = y0
            if x1 > X1: X1 = x1
            if y1 > Y1: Y1 = y1
                
            poly = Polygon(list(footprint), alpha=0.5, color=colors[i])
            ax1.add_patch(poly) 
        
        ax1.set_xlim(X0-0.1,X1+0.1)
        ax1.set_ylim(Y0-0.1,Y1+0.1)  
        
        plt.show()

    def plot_sources(self, x=None, y=None, labels=None, vmin=-0.05, vmax=0.05, mag_cut=17, 
                     plot_labels=True, plot_all=True):

        fig = plt.figure()    
        
        if plot_all:
            #print(labels.shape)
            filt = np.ones(labels.shape[0], dtype=int)
        else:
            filt = labels < mag_cut

        ax1 = fig.add_subplot(121)
        ax1.imshow(self.direct[0], origin="lower", vmin=vmin, vmax=vmax)
        ax1.scatter(x[filt], y[filt], fc="None", ec="r", s=60, lw=0.3)
        ax1.set_title("Direct")

        ax2 = fig.add_subplot(122)
        ax2.imshow(self.slitless[0], origin="lower", vmin=vmin, vmax=vmax)
        ax2.scatter(x[filt], y[filt], fc="None", ec="r", s=60, lw=0.3)
        ax2.set_title("Slitless")

        if plot_labels:
            for j,label in enumerate(labels):
                if label < mag_cut:
                    ax1.text(x[j], y[j], label)
                    ax2.text(x[j], y[j], label)
                    
        plt.show()

    def plot_sources_wcs(self, x=None, y=None, labels=None, vmin=-0.05, vmax=0.05, mag_cut=17, 
                     plot_labels=True, plot_all=True):

        fig = plt.figure()    
        
        if plot_all:
            #print(labels.shape)
            filt = np.ones(labels.shape[0], dtype=int)
        else:
            filt = labels < mag_cut

        ax1 = fig.add_subplot(121, projection=self.wcs_direct)
        ax1.imshow(self.direct[0], origin="lower", vmin=vmin, vmax=vmax)
        ax1.scatter(x[filt], y[filt], fc="None", ec="r", s=60, lw=0.3, 
                    transform=ax1.get_transform('world'))
        ax1.set_title("Direct")

        ax2 = fig.add_subplot(122, projection=self.wcs_slit)
        ax2.imshow(self.slitless[0], origin="lower", vmin=vmin, vmax=vmax)
        ax2.scatter(x[filt], y[filt], fc="None", ec="r", s=60, lw=0.3,
                    transform=ax2.get_transform('world'))
        ax2.set_title("Slitless")

        if plot_labels:
            for j,label in enumerate(labels):
                if label < mag_cut:
                    ax1.text(x[j], y[j], label, transform=ax1.get_transform('world'))
                    ax2.text(x[j], y[j], label, transform=ax2.get_transform('world'))
                    
        plt.show()

        
    def plot_sens(self):
        fig = plt.figure()

        ax1 = fig.add_subplot(111)
        ax1.plot(self.sens_cat['WAVELENGTH'], self.sens_cat['SENSITIVITY'],label="Roman Prism")
        ax1.set_xlim(6000,20000)
        ax1.set_ylabel("e-/s per erg/s/cm$^2$/Ang") # according to axe_manual
        ax1.legend()
        
        plt.show()


    def clip_bad_pixels(self, thresh=1):

        cutout = self.cutout[0]

    def spectral_trace(self, disp_axis=0, nseg=-1, vmin=-1, vmax=1, flip_axes=0):

        cutout = self.cutout[0]
        print(cutout.shape)

        N = cutout.shape[disp_axis]
        print(N)


        nl = np.arange(N)

        #nseg = 23
        #nseg = 56
        #nseg = 557
        if nseg == -1:
            nseg = N
        
        seg = int(N/nseg)
        print(seg)


        params = []
        param_errors = []


        nl_bin = []

        for i in range(nseg):
        #for i in range(N):

            if disp_axis == 0:
                #profile0 = cutout[i,:]
                profile0 = np.nansum(cutout[i*seg:(i+1)*seg, :], axis=0)
            
            elif disp_axis == 1:
                #profile0 = cutout[:,i]
                profile0 = np.nansum(cutout[:, i*seg:(i+1)*seg], axis=1)

            nl_bin.append(np.mean(nl[i*seg:(i+1)*seg]))

            p, perr = fit_gauss(profile0, verb=0, noise=0.05)

            params.append(p)
            param_errors.append(perr)

        nl_bin = np.array(nl_bin)

        params = np.array(params)
        param_errors = np.array(param_errors)

        print(params)
        print(param_errors)


        fig = plt.figure()

        ax1 = fig.add_subplot(211)

        if flip_axes==0:
            ax1.imshow(cutout, origin="lower", vmin=vmin, vmax=vmax)
        elif flip_axes==1:
            ax1.imshow(cutout.T, origin="lower", vmin=vmin, vmax=vmax)

        ax1.scatter(nl_bin, params[:,1], marker="x", c="tab:red")
        ax1.set_aspect(3.)

        ax2 = fig.add_subplot(212)
        ax2.errorbar(nl_bin, params[:,1], yerr=param_errors[:,1])

        ax2.set_xlabel("X pos (pixels)")
        ax2.set_ylabel("Y fit (pixels)")


        plt.show()
            
        fit_poly(nl_bin, params[:,1], param=[0, 1])

          
        
    def extract_cutout(self, disp_axis=0, spatial_axis=1, flip_axis=0):

        raw_sci = np.nansum(self.cutout[0], axis=disp_axis)
        raw_err = np.sqrt(np.nansum(self.cutout[1]**2, axis=disp_axis))

        print(self.wav[0], self.wav[-1])
        print(len(self.wav))
        print(raw_sci.shape)

        sens = np.array([self.interp_sens(w) for w in self.wav])

        if flip_axis:
            raw_sci = np.flip(raw_sci)
            raw_err = np.flip(raw_err)

        
        flux_sci = raw_sci/sens
        flux_err = raw_err/sens

        self.spectrum = [
            raw_sci,
            raw_err,
            flux_sci,
            flux_err,
        ]

        self.profile = np.nansum(self.cutout[0], axis=spatial_axis)

        return self.spectrum, self.profile


def zodi_subtract(sci_file, zodi_file, out_file, spec_exptime=574., spec_gain=2.):
    hdu_sci = pyfits.open(sci_file)
    sci = hdu_sci[1].data
    
    hdu = pyfits.open(zodi_file)
    zodi = hdu[1].data
    
    sci_sub = 1.0*sci - 1.0*zodi
        
    hdu_sci[1].data = sci_sub/spec_exptime/spec_gain
    
    chi2 = hdu_sci[2].data
    hdu_sci[2].data = np.sqrt(1.0*chi2)/spec_exptime/spec_gain
        
    hdu_sci.writeto(out_file, overwrite=True)


def insidePolygon(A,p):
    # format:

    #  points in polygon
    #  A = [[x0,y0],
    #       [x1,y1],
    #       [x2,y2]]

    #  test point
    #  p = [ x, y ]

    n = len(A)
    angle = 0
    for i in range(n):
        x1 = A[i][0] - p[0]
        y1 = A[i][1] - p[1]
        x2 = A[(i+1)%n][0] - p[0]
        y2 = A[(i+1)%n][1] - p[1]
        angle += angle2d(x1,y1,x2,y2)
    if np.abs(angle) < np.pi:
        return False
    else:
        return True

def angle2d(x1,y1,x2,y2):
    theta1 = np.arctan2(y1,x1)
    theta2 = np.arctan2(y2,x2)
    dtheta = theta2 - theta1
    while (dtheta > np.pi):
        dtheta -= 2*np.pi
    while (dtheta < -np.pi):
        dtheta += 2*np.pi
    return dtheta


def inspect_direct(ra0, dec0, cat, direct_file, plot_labels=1, hw_box=50, ext=1,
                   ra_key="RA", dec_key="DEC", mag_key="MAG_AUTO", vmin=-1, vmax=1):

    """not sure this works"""

    pf = pyfits.open(direct_file)
    head = pf[ext].header 
    img = pf[ext].data

    wcs = WCS(header=head)

    coords = [[ra0,dec0]]
    pixels = wcs.wcs_world2pix(coords, 1)
    #print(pixels)
    xc, yc = pixels[0]

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection=wcs)
    ax1.imshow(img, vmin=vmin, vmax=vmax, origin="lower")

    ax1.scatter(cat[ra_key], cat[dec_key], fc="None", ec="tab:red", s=50, lw=1.0, alpha=1,
                        transform=ax1.get_transform('world'))

    if plot_labels:
        for j,mag in enumerate(cat[mag_key]):
            ax1.text(cat[ra_key][j] - 0.3/3600., 
                     cat[dec_key][j] + 0.3/3600., 
                     "%.1f" % (mag), 
                     transform=ax1.get_transform('world'),
            )
    #print(ax1.get_xlim())
    #print(ax1.get_ylim())

    ax1.set_xlim(xc - hw_box, xc + hw_box)
    ax1.set_ylim(yc - hw_box, yc + hw_box)

    #print(ra0, dec0)
    
    plt.show()



def ab_to_fnu(mab=None):
    """Convert AB magnitude to Fnu [erg/s/cm^2/Hz]"""
    # mab = -2.5*log10(fnu) - 48.6
    fnu = 10 ** (-0.4 * (mab + 48.6))
    return fnu

def fnu_to_flam(fnu=None, wav=None):
    """Convert Fnu [erg/s/cm^2/Hz] to Flam [erg/s/cm^2/Ang]"""
    flam = fnu * c / (wav * Ang) ** 2 * Ang
    return flam

def flam_to_fnu(flam=None, wav=None):
    """Convert Flam [erg/s/cm^2/Ang to Fnu [erg/s/cm^2/Hz]]"""
    fnu = flam * (wav * Ang) ** 2 / c / Ang
    return fnu

#    # effective wavelength (pivot wavelength)
#    # http://www.astro.ljmu.ac.uk/~ikb/research/mags-fluxes/
#    numerator = integrate.trapz(tran,wav)
#    denominator = integrate.trapz(tran/wav**2,wav)
#    effwav_sq = numerator/denominator
#    effwav = np.sqrt(effwav_sq)

def calc_pivot(wav,tran):
    
    numerator = integrate.trapz(tran * wav, wav)
    denominator = integrate.trapz(tran / wav, wav)    
    pivot = np.sqrt(numerator/denominator)
    
    return pivot



def scale_spectra(spec=None, filt=None, mag=None, wav=None):
    """
    flux in erg/s/cm^2/Ang
    wavelength in Ang
    """
    
    spec_wav, spec_flux = spec
    filt_wav, filt_tran = filt
    
    w0, w1 = filt_wav[0], filt_wav[-1]
    i0 = np.argmin(np.abs(spec_wav - w0)) + 1
    i1 = np.argmin(np.abs(spec_wav - w1)) - 1
    #print(w0,w1)
    #print(i0,i1)

    interp_filt = interpolate.interp1d(filt_wav, filt_tran)
    filt_tran_rebin = np.array([interp_filt(w) for w in spec_wav[i0:i1]])
    filt_norm = np.sum(filt_tran_rebin)
    filt_max = np.max(filt_tran_rebin)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(spec_wav[i0:i1],filt_tran_rebin)
    plt.show()  
    
    filt_tran_rebin /= filt_norm
    
    #print(filt_norm)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(spec_wav[i0:i1],filt_tran_rebin)
    plt.show()  
        
    spec_flux_filt_conv = spec_flux[i0:i1]*filt_tran_rebin
    spec_flux_filt_norm = np.sum(spec_flux_filt_conv)
    
    #print(spec_flux_filt_norm)
    
    fnu = ab_to_fnu(mab=mag)   
    flam = fnu_to_flam(fnu=fnu, wav=wav)
    #print(flam)
    
    spec_flux_scale = spec_flux*flam/spec_flux_filt_norm*filt_norm
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.plot(filt_wav,filt_tran)
    #plt.show()   
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.plot(spec_wav[i0:i1], spec_flux[i0:i1])
    ax.plot(spec_wav[i0:i1], spec_flux_filt_conv)
    plt.show()   
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(spec_wav[i0:i1], spec_flux_scale[i0:i1])
    plt.show()
    
    #print(spec_flux_filt_norm)
    
    return spec_wav, spec_flux_scale


class EuclidData:
    def __init__(
        self,
        catalog1 = "TestPoints_highSNR_mod1_14324.fits",
        catalog2 = "NIS_detector_catalog_11.cat",
        model_spectra = "NIS_catalog_file_11.spc.fits",
        euclid_sims_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/Raw/EuclidSIMS/",
        catalog1_path = "./",
        catalog2_path = "TestPoints_highSNR_mod1_14324_2023_05_26/frame1/Catalogs/",
        model_spectra_path = "TestPoints_highSNR_mod1_14324_2023_05_26/frame1/Spectra/",
    ):
        
        catalog1_file = euclid_sims_path + catalog1_path + catalog1
        catalog2_file = euclid_sims_path + catalog2_path + catalog2
        
        self.model_spectra_file = euclid_sims_path + model_spectra_path + model_spectra
        
        ##########
        # catalog1
        ##########
        self.cat1 = Table.read(catalog1_file, format="fits")
        N1 = len(self.cat1)
        self.cat1["index"] = np.arange(N1)
        #print([col for col in self.cat1.colnames if "TU_" in col])
        
        self.Euclid_bands = ['VIS','NIR_Y','NIR_J','NIR_H']
        self.Euclid_bands_flux = ['TU_FNU_VIS_MAG', 'TU_FNU_Y_NISP_MAG', 'TU_FNU_J_NISP_MAG', 
                                  'TU_FNU_H_NISP_MAG'] 

        for bk,fk in zip(self.Euclid_bands, self.Euclid_bands_flux):
            fnu_Jy = self.cat1[fk] # Jy    
            mab = -2.5*np.log10(fnu_Jy) + 8.90   # Jy --> AB mag            
            mab[np.isinf(mab)]=-99.
            self.cat1[bk] = mab  
    
        #print(self.cat1.colnames)
        #print(len(self.cat1))

        ##########
        # catalog2
        ##########
        self.cat2 = Table.read(catalog2_file, format="ascii") 
        N2 = len(self.cat2)
        self.cat2["index"] = np.arange(N2)
        
        self.get_spectra_ids()

        
    def get_spectra_ids(self):
        t0 = time.time()

        with pyfits.open(self.model_spectra_file) as hdul:                 # 28 seconds
        #with pyfits.open(model_spectra, memmap=False) as hdul:
        #with pyfits.open(model_spectra, memmap=True) as hdul:
            N = len(hdul)
            #print(N)
            #print(hdul[0].header)
            #print(hdul[1].header)
            #print(hdul[1].columns)
    
            self.ids = [hdul[i].header['ID'] for i in range(1,N)]
        #print(ids)

        t1 = time.time()
        print("Time to compute: %.2f s" % (t1-t0))
        
        return self.ids
    
    def get_model_spectra(self, ra, dec):
        
        sid, _ = self.find_source(ra, dec)
        
        #print(sid) # SOURCE_ID
        ind = self.ids.index(sid)+1
        #print(ind) # index
        
        with pyfits.open(self.model_spectra_file) as hdul:
            print('index =', ind)
            head = hdul[ind].header
            data = hdul[ind].data
            self.model_wav = data.field('WAV_NM')*10 # Ang
            self.model_flux = data.field('FLUX') # ? 
            print('ID =', head['ID'])
            
        return self.model_wav, self.model_flux
    
    def plot_model_spectra(self):
        #print(self.model_wav)
        #print(self.model_flux)

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(self.model_wav, self.model_flux)
        plt.show()
        
    #def calc_euclid_mags_from_model(self):
    
    def find_source(self, ra, dec):
        dr1 = np.sqrt((self.cat1['RA_MAG']-ra)**2*np.cos(dec/180*np.pi)**2 + 
                      (self.cat1['DEC_MAG']-dec)**2)*3600.

        self.cat1['dr'] = dr1

        ind = np.argmin(dr1)
        sid = self.cat1['SOURCE_ID'][ind]
        #print(np.argmin(dr1))
        mab = self.cat1['NIR_H'][ind]
        print('SOURCE_ID:%d, mag=%.2f, dr=%.2f"' % (sid, mab, np.min(dr1)))
        return sid, ind


# In[11]:


def get_model_spectra_A(model_spectra, ind):

    with pyfits.open(model_spectra) as hdul:                 # 28 seconds
    #with pyfits.open(model_spectra, memmap=False) as hdul:
    #with pyfits.open(model_spectra, memmap=True) as hdul:
        #hdul.info()
        #N = len(hdul)
        #print(N)
        #print(hdul[0].header)
        #print(hdul[1].header)
        #print(hdul[2].header)
    
        model_wav = hdul[1].data
    
        cube = hdul[2].data
        #print(cube.shape)
    
        model_flux = cube[ind,:]

        return model_wav, model_flux
    
def plot_model_spectra_A(model_spectra, ind):
    
    with pyfits.open(model_spectra) as hdul:                 # 28 seconds
    #with pyfits.open(model_spectra, memmap=False) as hdul:
    #with pyfits.open(model_spectra, memmap=True) as hdul:
        #hdul.info()
        #N = len(hdul)
        #print(N)
        #print(hdul[0].header)
        #print(hdul[1].header)
        #print(hdul[2].header)
    
        model_wav = hdul[1].data
    
        cube = hdul[2].data
        #print(cube.shape)
    
        model_flux = cube[ind,:]

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(model_wav, model_flux)
        plt.show()


class EuclidFilters:
    
    def __init__(self):
        self.filter_path = "/Users/gwalth/data/Roman/grizli/tech_files/Euclid"
        self.filters = {}
        
        self.load_nisp_filters()
        self.load_vis_filter()
        
        # erg/s/cm^s/Hz  microJy, nJy
        self.Consts = [48.6, -23.9, -31.4]
        self.units  = ["Jy", "uJy", "nJy"]
        
    def log_fn(self, y, C):
        return -0.4*(y+C)
        
    def mag_fn(self, y, C):       
        return -2.5*np.log10(y) - C
    
        # F(Jy) = 10**23 * 10**(-AB+48.6)/2.5
        #              23-AB/2.5-48.6/2.5 = -0.4*(AB-8.9)
        #logf = -0.4*(y-8.9) # log10(F[Jy])

        # F(microJy) = 10**6 * 10**23 * 10**(-AB+48.6)/2.5
        #              6+23-AB/2.5-48.6/2.5 = -0.4*(AB-23.9)
        #logf = -0.4*(y-23.9) # log10(F[microJy])

        # F(nJy) = 10**9 * 10**23 * 10**(-(AB+48.6)/2.5)
        #          9+23-AB/2.5-48.6/2.5 = -0.4*(AB-31.4)
        #logf = -0.4*(y-31.4) # log10(F[nJy])

        #f = 10**(logf) # Jy
        
    def load_nisp_filters(self):
        L = glob.glob(self.filter_path + "/NISP-PHOTO-PASSBANDS-V1-*.dat")
        #print(L)
        
        for l in L:
            # NISP-PHOTO-PASSBANDS-V1-J_throughput.dat
            name = os.path.basename(l)
            band = name.split("_")[0][-1:]

            tbl = Table.read(l, format="ascii")
            #print(tbl)
            wav = tbl["WAVE"]
            trans = tbl["T_TOTAL"]
            
            self.filters['NISP_' + band] = {"wav": wav, "trans": trans}
    
    def load_vis_filter(self):
        #vis_file = "SpaceSegment.PLM.PLMTransmissionVISEOL.table"
        vis_file = "SpaceSegment.PLM.PLMTransmissionVISNominalEOL.table"
        det_file = "SpaceSegment.Instrument.VIS.MeanDetectorQuantumEfficiencyNominalBOL.table"
        
        vis_tbl = Table.read(os.path.join(self.filter_path,vis_file), format="ascii")
        
        # end-of-life
        EOL = vis_tbl['F1']
        degradation = 0.05 
        # beginning-of-life
        # EOL = (1. - degradation) * BOL
        BOL = EOL/(1. - degradation)

        interp_vis = interpolate.interp1d(vis_tbl['lambdaVIS'], BOL)

        det_tbl = Table.read(os.path.join(self.filter_path,det_file), format="ascii")
        #print(det_tbl)

        #wav = det_tbl['lambda']
        wav = np.arange(300, 1100, 2.)
        
        interp_qe = interpolate.interp1d(det_tbl['lambda'], det_tbl['QE'])
        trans = [interp_qe(w0)*interp_vis(w0) for w0 in wav]
        
        self.filters['VIS'] = {'wav': wav, 'trans': trans}
        
    def plot_filters(self):
        
        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        for filt in self.filters:
            wav = self.filters[filt]['wav']
            trans = self.filters[filt]['trans']

            ax1.plot(wav,trans)
            ax1.set_xlabel("Wavelength [nm]")
            ax1.set_ylabel("Throughput")
            
        plt.show()
        
    def calc_spectra_mag(self, band, fnu, wav, unit="Jy", area=0.785, plot=1):
        """area (default) = pi * r**2 (r=0.5 arcsec)"""
        ############
        # sed info #
        ############
        
        # interpolate source spectrum
        S = interpolate.interp1d(wav, fnu)
        
        filt_wav = self.filters[band]['wav']
        filt_trans = self.filters[band]['trans']
        
        if plot:
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(filt_wav, filt_trans, c="tab:orange")
            
            ax2 = ax1.twinx()
            tmp_flam = fnu_to_flam(fnu=fnu, wav=wav)

            ax2.plot(wav, tmp_flam, c="tab:blue")
            
            plt.show()
        
        # interpolate filter
        F = interpolate.interp1d(filt_wav, filt_trans)

        filt_all_int = integrate.trapz(filt_trans, filt_wav)
        dw = wav[1]-wav[0]
        #filt_all_sum  = np.sum(tran)*dw
        #print(filt_all_int)
        #print filt_all_sum

        # true_flux = flux_measured_in_filter * C 
        F_all = np.array([S(w)*F(w) for w in filt_wav])

        flux_int = integrate.trapz(F_all, filt_wav)
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
        
        ind = self.units.index(unit)
        Const = self.Consts[ind]

        mAB = self.mag_fn(flux_corr, Const)
        print("mAB =", mAB)
        mAB_5sig = self.mag_fn(flux_corr/5., Const)
        print("mAB (5 sig) =", mAB_5sig)
        mAB_10sig = self.mag_fn(flux_corr/10., Const)
        print("mAB (10 sig) =", mAB_10sig)
        mAB_area = self.mag_fn(flux_corr/area, Const)
        print('mAB (mag/") =', mAB_area)
        print()


verb = 0

if verb:
    # https://www.hamamatsu.com/us/en/resources/interactive-tools/photon-flux-to-radiant-flux-calculator.html
    
    ####################################
    # constants
    ####################################
    Jy = 1.0E-23        # erg/s/cm^2/Hz
    Ang = 1E-8          # cm
    nm = 1E-7           # cm
    c = 2.9979E10       # cm/s
    #h = 6.626068E-27    # erg * s 
    h = 6.626068E-34    # J * s
    
    ####################################
    wav = 500 # nm
    E_phot = (h*c)/(wav*nm) # J
    print(E_phot)
    
    flux_phot = 1. # photons/s
    flux  = E_phot * flux_phot # J/s
    print(flux)
    
    
    # In[ ]:
    
    
    d = 1.2 # m
    collect_area = np.pi*(d*100./2)**2 # cm**2
    print(collect_area)
    
    
    # In[ ]:
    
    
    d = 7. # mm
    collect_area = np.pi*((d/10.)/2)**2 # cm**2
    print(collect_area)
    
    
    # In[ ]:
    
    
    c = 2.9979E10       # cm/s
    h = 6.626068E-27    # erg * s 
    
    #ABmag = 16.
    ABmag = 8.5
    #wav_cen = 17714. # Ang
    wav_cen = 5600 # Ang
    #wav_width = 4999. # Ang
    wav_width = 1000
    #################################################################
    fnu = 10**(-0.4*(ABmag + 48.60))                 # erg/s/cm^2/Hz
    print("Calculated from magnitude")
    print(fnu,"erg/s/cm^2/Hz")
    #flambda = fnu*Ang/((lam_obs*mu)**2/c)
    flambda = fnu*Ang/((wav_cen*Ang)**2/c)
    print(flambda,"erg/s/cm^2/Ang")
    print(flambda*wav_width,"erg/s/cm^2")
    print(flambda*wav_width*collect_area,"erg/s")
    E_phot = (h*c)/(wav_cen*Ang) # erg
    print(E_phot,"erg")
    print(flambda*wav_width/E_phot,"photons/s/cm^2")
    print(flambda*wav_width/E_phot*collect_area,"photons/s")
    print(flambda/E_phot,"photons/s/cm^2/Ang")
    print(flambda/E_phot*collect_area,"photons/s/Ang")
    # above is correct!!
    ########################################################################
    #photlam # photons/s/cm^2/Ang
    
    
    # ### Spectrum conversion test
    
    # In[ ]:
    
    
    datafiles = "/Users/gwalth/python.linux/datafiles"
    #spectrum = "vega_all.fits"
    spectrum = "spec_vega.fits"
    ext = 0
    spec_file = os.path.join(datafiles, spectrum)
    pf = pyfits.open(spec_file)
    
    # Vega test spectrum
    if spectrum == "vega_all.fits":
    
        spec = pf[ext].data
        head = pf[ext].header
        cdelt1 = head["cdelt1"]
        crval1 = head["crval1"]
    
        nelem = spec.shape[0]
        specwave = (np.arange(nelem))*cdelt1 + crval1  # Angstrom
        #spec /= 206265.**2
    
    elif spectrum == "spec_vega.fits":
    
        specwave = pf[ext].data[0,:] # Angstrom
        spec = pf[ext].data[1,:]     # erg/s/cm^2/Ang
        nelem = spec.shape[0]
    
    print(specwave)
    print(spec)
    print(nelem)
    
    E_phot = (h*c)/(specwave*Ang) # erg
    #print E_phot
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(specwave, spec/E_phot) # photons/s/cm^2/Ang
    ax.set_xlim(3000,11000)
    #p.set_xlim(20000,24000)
    ax.set_ylim(0,2000)
    ax.set_xlabel("Wavelength ($\AA$)")
    ax.set_ylabel("Flux (photons cm$^{-2}$ s$^{-1}$ $\AA^{-1}$)")
    ax.set_title("Vega photon spectrum")
    
    # ABnu = STlamb @ 5492.9 Ang
    STlamb = 3.63E-9*np.ones((nelem)) # erg/cm^2/s/Ang   STlamb = 0
    ABnu = 3.63E-20*np.ones((nelem))  # erg/cm^2/s/Hz    ABnu = 0
    
    ax.plot(specwave,STlamb/E_phot)
    ax.plot(specwave,ABnu/E_phot*(c/(specwave*Ang)**2)*Ang)
    
    
    plt.show()
    
    
    # ![STmag_ABmag_synphot.png](attachment:STmag_ABmag_synphot.png)
    
    # In[ ]:
    
    
    #PHOTMODE= 'WFC3 IR F140W'      / observation con                                
    #PHOTFLAM=        1.4737148E-20 / inverse sensitivity, ergs/cm2/Ang/electron     
    #PHOTFNU =        9.5291135E-08 / inverse sensitivity, Jy*sec/electron           
    #PHOTZPT =       -2.1100000E+01 / ST magnitude zero point                        
    #PHOTPLAM=        1.3922907E+04 / Pivot wavelength (Angstroms)                   
    #PHOTBW  =        1.1323900E+03 / RMS bandwidth of filter plus detector 
    
    
    # In[ ]:
    
    
    
    # ## Transmission
    
    # ### VIS
    
    # In[ ]:
    
    
    print(os.getcwd())
    
    from scipy import interpolate
    
    filter_path = "/Users/gwalth/data/Roman/grizli/tech_files/Euclid"
    L = glob.glob(filter_path + "/*.table")
    print(L)
    
    #vis_f = "SpaceSegment.PLM.PLMTransmissionVISEOL.table"
    vis_f = "SpaceSegment.PLM.PLMTransmissionVISNominalEOL.table"
    det_f = "SpaceSegment.Instrument.VIS.MeanDetectorQuantumEfficiencyNominalBOL.table"
    
    vis_tbl = Table.read(os.path.join(filter_path,vis_f), format="ascii")
    print(vis_tbl)
    #interp_vis = interpolate.interp1d(vis_tbl['lambdaVIS'],vis_tbl['Transmission'])
    interp_vis = interpolate.interp1d(vis_tbl['lambdaVIS'],vis_tbl['F1'])
    
    det_tbl = Table.read(os.path.join(filter_path,det_f), format="ascii")
    print(det_tbl)
    
    #wav = det_tbl['lambda']
    
    wav = np.arange(300,1100,2.)
    
    interp_qe = interpolate.interp1d(det_tbl['lambda'],det_tbl['QE'])
    
    trans = [interp_qe(w0)*interp_vis(w0) for w0 in wav]
    
    t_max = np.max(trans)
    ind = np.argmax(trans)
    
    print(t_max/2.)
    print(t_max)
    print(wav[ind])
    
    w0 = 560
    w1 = 900
    
    i0 = np.argmin(np.abs(wav-w0))
    i1 = np.argmin(np.abs(wav-w1))
    
    print(np.mean(trans[i0:i1]))
    print(np.median(trans[i0:i1]))
    print(np.std(trans[i0:i1]))
    
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(wav,trans)
    ax1.plot([300,1100],[t_max/2.,t_max/2.])
    #ax1.plot(vis_tbl['lambdaVIS'],vis_tbl['Transmission'])
    #ax1.plot(vis_tbl['lambdaVIS'],vis_tbl['F1'])
    #ax1.plot(det_tbl['lambda'],det_tbl['QE'])
    ax1.set_xlabel("Wavelength [nm]")
    ax1.set_ylabel("Throughput")
    plt.show()
    
    
    # ### NISP
    
    # In[ ]:
    
    
    np.pi*(0.5**2)
    
    
    # In[ ]:
    
    
    ef = EuclidFilters()
    ef.plot_filters()
    print(ef.filters.keys())
    
    
    # In[ ]:
    
    
    print(os.getcwd())
    
    # https://www.euclid.caltech.edu/page/filters
    
    # M. Schirmer et al. 2022
    # Euclid preparation
    # XVIII. The NISP photometric system
    
    # https://euclid.esac.esa.int/msp/refdata/nisp/NISP-PHOTO-PASSBANDS-V1
    # T_{\textrm{Tot}} = T_{\textrm{Tel}} * T_{\textrm{NIOA}} * T_{\textrm{Filt}} * T_{\textrm{QE}}
    
    filter_path = "/Users/gwalth/data/Roman/grizli/tech_files/Euclid"
    L = glob.glob(filter_path + "/*.dat")
    print(L)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    for l in L:
        tbl = Table.read(l, format="ascii")
        #print(tbl)
        wav = tbl["WAVE"]
        tran = tbl["T_TOTAL"]
        
        pivot = calc_pivot(wav,tran)
        print("Pivot =", pivot)
        
        ax1.plot(wav,tran)
        ax1.set_xlabel("Wavelength [nm]")
        ax1.set_ylabel("Throughput")
    
         
    plt.show()
    
    
    # ### VIS + NISP
    
    # In[ ]:
    
    
    print(os.getcwd())
    
    from scipy import interpolate
    
    filter_path = "/Users/gwalth/data/Roman/grizli/tech_files/Euclid"
    L = glob.glob(filter_path + "/*.table")
    print(L)
    
    #vis_f = "SpaceSegment.PLM.PLMTransmissionVISEOL.table"
    vis_f = "SpaceSegment.PLM.PLMTransmissionVISNominalEOL.table"
    det_f = "SpaceSegment.Instrument.VIS.MeanDetectorQuantumEfficiencyNominalBOL.table"
    
    vis_tbl = Table.read(os.path.join(filter_path,vis_f), format="ascii")
    #print(vis_tbl)
    #interp_vis = interpolate.interp1d(vis_tbl['lambdaVIS'],vis_tbl['Transmission'])
    
    # end-of-life
    EOL = vis_tbl['F1']
    degradation = 0.05 
    # beginning-of-life
    # EOL = (1. - degradation) * BOL
    BOL = EOL/(1. - degradation)
    
    interp_vis = interpolate.interp1d(vis_tbl['lambdaVIS'],BOL)
    
    
    det_tbl = Table.read(os.path.join(filter_path,det_f), format="ascii")
    #print(det_tbl)
    
    #wav = det_tbl['lambda']
    
    wav = np.arange(300,1100,2.)
    
    interp_qe = interpolate.interp1d(det_tbl['lambda'],det_tbl['QE'])
    
    trans = [interp_qe(w0)*interp_vis(w0) for w0 in wav]
    
    
    
    filter_path = "/Users/gwalth/data/Roman/grizli/tech_files/Euclid"
    L = glob.glob(filter_path + "/*.dat")
    print(L)
    
    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(111)
    
    ax1.plot(wav,trans)
    pivot = calc_pivot(wav,trans)
    print("Pivot =", pivot)
    
    for l in L:
        tbl = Table.read(l, format="ascii")
        print(tbl.columns)
        wav = tbl["WAVE"]
        tran = tbl["T_TOTAL"]
        
        pivot = calc_pivot(wav,tran)
        print("Pivot =", pivot)
        
        #ax1.plot(wav,tran)
        ax1.plot(wav,tbl["T_TOTAL"])
        #ax1.plot(wav,tbl["T_FILTER"])
        #ax1.plot(wav,tbl["T_TELESCOPE"])
        #ax1.plot(wav,tbl["T_NI-OA"])
        #ax1.plot(wav,tbl["T_QE"])
    
        ax1.set_xlabel("Wavelength [nm]")
        ax1.set_ylabel("Throughput")
    
         
    plt.show()
    
    
    # ### Calculate Pivot
    
    # In[ ]:
    
    
    print(os.getcwd())
    
    filter_path = "/Users/gwalth/data/Roman/grizli/tech_files/Euclid"
    L = glob.glob(filter_path + "/*.dat")
    print(L)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    
    for l in L:
        tbl = Table.read(l, format="ascii")
        #print(tbl)
        wav = tbl["WAVE"]
        tran = tbl["T_TOTAL"]
        
        pivot = calc_pivot(wav,tran)
        print("Pivot =", pivot)
        
        ax1.plot(wav,tran)
        ax1.set_xlabel("Wavelength [nm]")
        ax1.set_ylabel("Throughput")
    
         
    plt.show()
    
    
    # In[ ]:
    
    
    
    
