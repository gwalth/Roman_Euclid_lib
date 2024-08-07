import sys
import matplotlib.pyplot as plt
from grizli_functions import display_grizli
from grizli_aXe_glw import EuclidData

# usage:
#    python test_display1d.py FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11 DET11 156



root = sys.argv[1]
det_id = sys.argv[2]
id = int(sys.argv[3])
frame = 1

det = "DET" + det_id

print(len(sys.argv))
if len(sys.argv) > 4:
  z_in = float(sys.argv[4])
else:
  z_in = 0.0


euc_dat = EuclidData(
    catalog1 = "FSpatch_mod3_16183_TAcalib_newGRID_V1.fits",
    catalog2 = "NisInputConfiguration_%i_catalog_%s.cat" % (frame, det_id),
    model_spectra = "catalog_%s_%s.%i.spc.fits" % (det_id, det_id, frame),
    euclud_sims_path = "/Users/gwalth/data/Roman/grizli/sims/Euclid/Raw/EuclidSIMS/",
    catalog1_path = "FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/",
    catalog2_path = "FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/frame%i/Catalogs/" % (frame),
    model_spectra_path = "FSpatch_mod3_16183_TAcalib_V1_R_newGRID_2024-05-17/frame%i/Input_Spectra/" % (frame),
)





vmin = None
vmax = None

# FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11_DET11_00156
prefix = '{0}_{1}_{2:05d}'.format(root, det, id)

fig = plt.figure()
ax, info = display_grizli(prefix, w0=1.15, w1=1.95, labels=1, y0=-1, y1=-1, z_in=z_in, path="", lw=2, 
    fontsize=8, dispersers=['RED'], yscale=1.2, fig=fig, cmap='viridis_r',
    scale_size=1, vmin=vmin, vmax=vmax)


###########################################################################
sid, ind = euc_dat.find_source(info["RA"], info["DEC"])
model_wav, model_flux = euc_dat.get_model_spectra(info["RA"], info["DEC"])

w0, w1 = ax[3].get_xlim()
y0, y1 = ax[3].get_ylim()

ax[3].plot(model_wav/1e4, model_flux, c="tab:orange", label="input")

ax[3].set_xlim(w0, w1); ax[3].set_ylim(y0,y1)
ax[3].legend()
###########################################################################


plt.show()
