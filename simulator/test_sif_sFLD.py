import sys
sys.path.append("../")

import numpy as np
from matplotlib import pyplot as plt
from ibis_simulator import libRadSpectra

from proc_ibis import sif_sFLD, sif_3FLD

if __name__=="__main__":

    ibis_wavls=np.genfromtxt("ibis_wavelengths.txt")    


    w=libRadSpectra("white_ref_z1000m.out",ftype="TXT",dataCol=2,hdrLines=0)
    w.convert_units_Qcm2_to_Wm2()
    w.resample_to_ibis(ibis_wavls)

    f=libRadSpectra("uvspec_fluorescence_z1000m.out" ,ftype="TXT",dataCol=2,hdrLines=0)
    f.resample_to_ibis(ibis_wavls)
    f.convert_units_Qcm2_to_Wm2()

    f_wvl_o2a=f.wavl_of_min_over_band((759,770))

    a=libRadSpectra("../../libRadtran-2.0.4/examples/UVSPEC_FLUORESCENCE.FLU" ,ftype="TXT",dataCol=1,hdrLines=6)
    a.convert_units_Qcm2_to_Wm2()
    idx=np.abs(a.wavl-f_wvl_o2a).argmin()


    print(sif_sFLD(f,w)*1000)
    print(sif_3FLD(f,w)*1000)
    print(a.data[idx]*1000)
        
    #plt.plot(w.wavl,f.data/w.data)
    #plt.show()
    
    

