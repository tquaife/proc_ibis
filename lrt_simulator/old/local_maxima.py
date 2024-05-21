import sys
sys.path.append("../")

import numpy as np
import matplotlib.pyplot as plt
from ibis_simulator import libRadSpectra, convert_units_Qcm2_to_Wm2
from sif_retrieval import *

fn="lrt_io/uvspec--sza_0p0--umu_1p0--phi_0p0--zout_1p0--fluorescence_fluorescence_0p0--albedo_albedo_0p1.out"

ibis_wavls=np.genfromtxt("ibis_wavelengths.txt")   

l=libRadSpectra(fn)
l.resample_to_ibis(ibis_wavls)
l.convert_units_Qcm2_to_Wm2()

wl_beg=740
wl_end=745
m=get_local_maxima(l,wl_beg,wl_end)
n=get_local_minima(l,wl_beg,wl_end)

plt.xlim([wl_beg,wl_end])
plt.plot(l.wavl,l.data)
plt.plot(m.wavl,m.data,"o")
plt.plot(n.wavl,n.data,"o")
plt.show()


