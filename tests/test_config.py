##############################################################
# model.py

# line 4160

# conf_args = dict(instrume=self.grism.instrument,
#                  filter=direct_filter,
#                  grism=self.grism.filter,
#                  module=self.grism.module,
#                  chip=self.grism.ccdchip)

# self.conf_file = grismconf.get_config_filename(**conf_args)
# conf = grismconf.load_grism_config(self.conf_file)
##############################################################

from grizli import grismconf

conf_args = dict(instrume = "NIRISS",
                 grism = "GR150R",
                 filter = "F115W")


conf_file = grismconf.get_config_filename(**conf_args)
print(conf_file)
conf = grismconf.load_grism_config(conf_file)
print(conf)
print(dir(conf))
print(conf.__dict__)
print(conf.conf.keys())


print()
print()
print()

sys.exit()


conf_args = dict(instrume = "NISP",
                 grism = "RGS000p0",
                 chip = "11")


conf_file = grismconf.get_config_filename(**conf_args)
print(conf_file)
conf = grismconf.load_grism_config(conf_file)
print(conf)
print(dir(conf))
print(conf.__dict__)
print(conf.conf.keys())

#conf_file = "../grizli/CONF/Euclid/GLWv3/RGS000_0/config_axesim_det_11.conf"
#conf_file = "../grizli/CONF/G141.F140W.V4.32.conf"
#conf = grismconf.load_grism_config(conf_file)
#print(conf)
#print(dir(conf))
#print(conf.__dict__)
#print(conf.conf.keys())

#for k in ['INSTRUMENT', 'CAMERA', 'GFILTER', 'DFILTER']:
#    print(conf.conf[k])


#print()
#print()
#print()
#print()
#print()

#instrume = 'NISP-GLWv3-det11'
#print(instrume.split("-"))

#print()
#print()
#print()
#print()
#print()
#det = instrume.split("-")[-1][-2:]
#print(det)

#print()
#print()
#print()
#print()
##print()
#conf_file = grismconf.get_config_filename(instrume='NISP', grism='RGS000p0', chip='11', module=None)
#print(conf_file)
#conf = grismconf.load_grism_config(conf_file)
#print(conf)
#print(dir(conf))
#print(conf.__dict__)


#print()
#print()
#print()
#print()
#print()



#conf_file = grismconf.get_config_filename(instrume='WFC3', filter='F140W', grism='G141', module=None, chip=1)
#print(conf_file)
#conf = grismconf.load_grism_config(conf_file)
#print(conf)
#print(dir(conf))
#print(conf.__dict__)

#print()
#print()
#print()
#print()
#print()


#conf_file = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Roman.G150-v2-GLW.conf"
#conf = grismconf.load_grism_config(conf_file)
#print(conf)
#print(dir(conf))
#print(conf.__dict__)

#print()
#print()
#print()
#print()
#print()

#conf_file = "/Users/gwalth/data/Roman/grizli/grizli/CONF/Euclid/GLWv3/RGS000p0/config_axesim_det_11_test.conf"
#conf = grismconf.load_grism_config(conf_file)
#print(conf)
#print(dir(conf))
#print(conf.__dict__)

#print()
#print()
#print()
#print()
#print()
#conf_file = "/Users/gwalth/data/Roman/grizli/grizli/CONF/GR150C.F090W.conf"
#conf = grismconf.load_grism_config(conf_file)
#print(conf)
#print(dir(conf))
#print(conf.__dict__)
