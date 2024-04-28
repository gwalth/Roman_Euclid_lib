import sys
import grizli
import numpy
import scipy
import astropy

print('\n Python version: ', sys.version)
print('\n Numpy version: ', numpy.__version__)
print(numpy.__file__)
print('\n Scipy version: ', scipy.__version__)
print(scipy.__file__)
print('\n Grizli version: ', grizli.__version__)
print(grizli.__file__)
print('\n Astropy version: ', astropy.__version__)
print(astropy.__file__)



#import importlib
#import sys

## Grizli and requirements
#import grizli

### Module versions
#print(sys.version + '\n')

#for module in ['grizli','grizli_aws', 'eazy', 'reprocess_wfc3', 'tristars', 'mastquery', 
#               'wfc3dash', 'prospect', 'sep', 'numpy', 'scipy', 'astropy', 'astroquery', 
#               'shapely', 'photutils', 'drizzlepac', 'wfc3tools', 'stsci.tools']:
#    #print(module)
#    try:
#        mod = importlib.import_module(module)
#        print('{0:>20} : {1}'.format(module, mod.__version__))
#    except ModuleNotFoundError:
#        print('{0:>20} : {1}'.format(module, '*failed*'))
