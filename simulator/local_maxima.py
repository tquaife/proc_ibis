import numpy as np
import matplotlib.pyplot as plt
from ibis_simulator import libRadSpectra, convert_units_Qcm2_to_Wm2

fn="lrt_io/uvspec--sza_0p0--umu_1p0--phi_0p0--zout_1p0--fluorescence_fluorescence_0p0--albedo_albedo_0p1.out"

ibis_wavls=np.genfromtxt("ibis_wavelengths.txt")   

l=libRadSpectra(fn)
l.resample_to_ibis(ibis_wavls)
l.convert_units_Qcm2_to_Wm2()

wl_beg=740
wl_end=745
idx_beg = int(np.abs(l.wavl-wl_beg).argmin())
idx_end = int(np.abs(l.wavl-wl_end).argmin())

m=libRadSpectra()
for i in range(int(idx_beg),int(idx_end)):
    if l.data[i]>l.data[i-1] and l.data[i]>l.data[i+1]:
        m.data=np.append(m.data,l.data[i])
        m.wavl=np.append(m.wavl,l.wavl[i])

plt.xlim([wl_beg,wl_end])
plt.plot(l.wavl,l.data)
plt.plot(m.wavl,m.data,"o")
plt.show()


