import numpy as np
from spectra_tools.spectra import Spectra
from sif_retrieval import *

import matplotlib.pyplot as plt

if __name__=="__main__":

    wref=Spectra(fname="./data_in/white_reference.csv",ftype="CSV",hdrLines=0)
    test=Spectra(fname="./data_in/test_spectra.csv",ftype="CSV",hdrLines=0)

    sif_opts_FLD_Cendrero19_O2b.in_wband=(686, 688) 
    sif_opts=sif_opts_FLD_Cendrero19_O2b
      
    print(sif_sFLD(test,wref,sif_opts=sif_opts)*0.01)
    print(sif_3FLD(test,wref,sif_opts=sif_opts)*0.01)
    print(sif_linReg(test,wref,sif_opts=sif_opts)*0.01)

    pseudo=sif_linReg_compute_pseudo_inverse(wref, sif_opts)
    print(sif_linReg_use_precomp_inverse(test, pseudo, sif_opts)*0.01)


    #plt.plot(wref.wavl,wref.data)
    #plt.plot(test.wavl,test.data)
    #plt.show()

