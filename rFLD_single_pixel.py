from proc_ibis import ibisBilReader
from sif_retrieval import *

import matplotlib.pyplot as plt

if __name__=="__main__":

    line=2600
    sample=100
    
    #read in reference data
    wref=Spectra(fname="./data_in/white_reference.csv",ftype="CSV",hdrLines=0)

    #read ibis data
    ibis_filename="./data_in/i200073b_mapped_utm30n.bil"
    ibis_data=ibisBilReader(ibis_filename)    
    ibis_spectra=ibis_data.get_pixel(sample,line)
    
    sif_opts=sif_opts_FLD_Cendrero19_O2a
    sif_opts.in_wband=(755, 770)
        
    gR=100000.0
    gF=500000.0
    oR=3
    oF=2
    
    pseudo=sif_ridgeReg_compute_pseudo_inverse_RFopts(wref,orderR=oR,gammaR=gR,orderF=oF,gammaF=gF,sif_opts=sif_opts)
    K=sif_ridgeReg_compute_pseudo_inverse_RFopts(wref,orderR=oR,gammaR=gR,orderF=oF,gammaF=gF,sif_opts=sif_opts,return_k=True)
    #F=sif_ridgeReg_use_precomp_inverse(ibis_spectra,pseudo,sif_opts=sif_opts, out_wvl=761)*ibis_data.units_to_mw_m2
    (R,F)=sif_ridgeReg_use_precomp_inverse_full_output(ibis_spectra,pseudo,sif_opts=sif_opts)

    F.data*=ibis_data.units_to_mw_m2
    print(F.closest_to_wavl(761)[1])
    
    ibis_spectra.trim(sif_opts.in_wband[0],sif_opts.in_wband[1])
    plt.plot(ibis_spectra.wavl,ibis_spectra.data)

    plt.twinx()
    plt.plot(F.wavl,F.data,"g-")
    plt.plot(R.wavl,R.data,"r-")
    plt.xlabel("wavelength")
    plt.show()
