
import glob, sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

# usage:
#    python ../../../../notebooks/test_inspect2d.py "*stack.fits"
#    python ../../../../notebooks/test_inspect2d.py "*beams.fits"


def arr_stats(arr, verb=1):

    X = arr.flatten()

    x0 = np.nanmin(X)
    x1 = np.nanmax(X)
    x_mean = np.nanmean(X)
    x_median = np.nanmedian(X)
    x_std = np.nanstd(X)
    q_2_3 = np.quantile(X, 0.023)
    q_15_9 = np.quantile(X, 0.159)
    q_84_1 = np.quantile(X, 0.841)
    q_97_7 = np.quantile(X, 0.977)

    if verb:
        print("min = %f" % (x0))
        print("max = %f" % (x1))
        print("mean = %f" % (x_mean))
        print("median = %f" % (x_median))
        print("std = %f" % (x_std))
        print("q(2.3) = %f" % (q_2_3))
        print("q(15.9) = %f" % (q_15_9))
        print("q(84.1) = %f" % (q_84_1))
        print("q(97.7) = %f" % (q_97_7))

    all_dict = {
        "min": x0,
        "max": x1,
        "mean": x_mean,
        "median": x_median,
        "std": x_std,
        "2.3": q_2_3,
        "15.9": q_15_9,
        "84.1": q_84_1,
        "97.7": q_97_7,
    }

    return all_dict
    



f_str = sys.argv[1]
print(f_str)

L = glob.glob(f_str)
L.sort()
print(L)

N = len(L)

cols = 2
rows = int(np.ceil(N/cols))

print(cols)
print(rows)




#fig, axs = plt.subplots(rows, cols, figsize=(16, 1.5), layout='constrained')
fig, axs = plt.subplots(rows, cols, figsize=(4*cols, 0.375*rows), layout='constrained')

for ax,l in zip(axs.flat,L):

    print(l)

    pf = pyfits.open(l)
    #pf.info()

    img = pf['SCI'].data
    head = pf['SCI'].header

    det_id = head["DET_ID"]

    stat_dict = arr_stats(img)

    std = stat_dict["std"]
    med = stat_dict["median"]
    vmin = med - 0.1*std
    vmax = med + 0.1*std

    ax.imshow(img, vmin=vmin, vmax=vmax)
    ax.text(0.025, 0.5, "DET" + det_id, transform=ax.transAxes)
    ax.set_axis_off()


    #ax.set_title(f'markevery={markevery}')
    #ax.plot(x, y, 'o', ls='-', ms=4, markevery=markevery)

plt.show()
