import sys
import matplotlib.pyplot as plt
from grizli_functions import display_grizli

# usage:
#    python test_display1d.py FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11 DET11 156


root = sys.argv[1]
det = sys.argv[2]
id = int(sys.argv[3])

print(len(sys.argv))
if len(sys.argv) > 4:
  z_in = float(sys.argv[4])
else:
  z_in = 0.0

# FSpatch_mod3_16183_TAcalib_V1_RNDC_2024-03-11_DET11_00156
prefix = '{0}_{1}_{2:05d}'.format(root, det, id)

fig = plt.figure()
ax = display_grizli(prefix, w0=1.15, w1=1.95, labels=1, y0=-1, y1=-1, z_in=z_in, path="", lw=2, 
                   fontsize=8, dispersers=['RED'], yscale=1.2, fig=fig, cmap='viridis_r',
                   scale_size=1)

plt.show()
