from grizli import grismconf
import matplotlib.pyplot as plt

#conf_file = "../grizli/CONF/Euclid/GLWv1/NISP_RGS000_21.conf"
#conf = grismconf.aXeConf(conf_file)
##conf.get_beams()
#fig = conf.show_beams()

#plt.show()


conf_file = "../../grizli/CONF/Euclid/GLWv3/RGS000p0/config_axesim_det_11.conf"
conf = grismconf.aXeConf(conf_file)
##conf.get_beams()
fig1 = conf.show_beams()

conf_file = "../../grizli/CONF/Euclid/GLWv3/RGS000m4/config_axesim_det_11.conf"
conf = grismconf.aXeConf(conf_file)
#conf.get_beams()
fig2 = conf.show_beams()

conf_file = "../../grizli/CONF/Euclid/GLWv3/RGS180p0/config_axesim_det_11.conf"
conf = grismconf.aXeConf(conf_file)
#conf.get_beams()
fig3 = conf.show_beams()

conf_file = "../../grizli/CONF/Euclid/GLWv3/RGS180p4/config_axesim_det_11.conf"
conf = grismconf.aXeConf(conf_file)
#conf.get_beams()
fig4 = conf.show_beams()


plt.show()
#print(conf)
