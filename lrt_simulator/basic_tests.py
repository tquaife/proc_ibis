import sys
sys.path.append("../")

from copy import copy, deepcopy

import numpy as np
from matplotlib import pyplot as plt

from ibis_simulator import libRadSpectra, convert_units_Qcm2_to_Wm2
from sif_retrieval import sif_sFLD, sif_3FLD, sif_opts_FLD_Cendrero19_O2a, sif_opts_FLD_Cendrero19_O2b
from run_libradtran import uvspec

if __name__=="__main__":

    #load ibis wavelengths
    ibis_wavls=np.genfromtxt("ibis_wavelengths.txt")    

    dataCol=1
    overwriteLrtOuts=False

    #generate reference spectra:
    u=uvspec()
    u.options["fluorescence"]="fluorescence 0.0"
    u.options["albedo"]="albedo 1.0"
    u.options["zout"]="sur"
    white_ref=libRadSpectra(u.run(overwrite=overwriteLrtOuts),dataCol=dataCol)
    white_ref.resample_to_ibis(ibis_wavls)
    white_ref.convert_units_Qcm2_to_Wm2()

    #create some SIF data with a spectrally flat
    #albedo and fluorescence
    F=0.015
    alb=0.1
    flat_sif=deepcopy(white_ref)
    flat_sif.data*=alb
    flat_sif.data+=F    

    #check the inversion algorithms
    #(they should all solve this problem perfectly)
    print("F: %0.3f alb: %0.3f"%(F*1000,alb))
    print("sFLD O2a:",sif_sFLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("sFLD O2b:",sif_sFLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)
    print("3FLD O2a:",sif_3FLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("3FLD O2b:",sif_3FLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)

    #plt.plot(white_ref.wavl,white_ref.data,label="alb=1.0")
    #plt.plot(flat_sif.wavl,flat_sif.data,label="alb=0.1")
    #plt.plot(flat_sif.wavl,flat_sif.data/white_ref.data,label="refl")
    #plt.legend()
    #plt.show()




