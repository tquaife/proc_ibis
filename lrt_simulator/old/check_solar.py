import numpy as np
from matplotlib import pyplot as plt
from ibis_simulator import libRadSpectra

if __name__=="__main__":

    a=libRadSpectra("solar_spectrum.txt",ftype="TXT",dataCol=1,hdrLines=0)
    a.convert_units_Qcm2_to_Wm2()
    b=libRadSpectra("../../libRadtran-2.0.4/data/solar_flux/atlas_plus_modtran",ftype="TXT",dataCol=1,hdrLines=0)
 
    w=libRadSpectra("lrt_io/white_ref_z1000m.out",ftype="TXT",dataCol=2,hdrLines=0)
    w.convert_units_Qcm2_to_Wm2()
    f=libRadSpectra("lrt_io/uvspec_fluorescence_z1000m.out" ,ftype="TXT",dataCol=2,hdrLines=0)
    f.convert_units_Qcm2_to_Wm2()

    
    plt.plot(a.wavl,a.data*1000)
    plt.plot(b.wavl,b.data)
    plt.plot(w.wavl,w.data*1000)
    plt.plot(f.wavl,f.data*1000)

    plt.show()
