import sys
sys.path.append("../")

import numpy as np
from matplotlib import pyplot as plt

from ibis_simulator import libRadSpectra, convert_units_Qcm2_to_Wm2
from sif_retrieval import sif_sFLD, sif_3FLD, sif_opts_FLD_Cendrero19_O2a, sif_opts_FLD_Cendrero19_O2b
from run_libradtran import uvspec


if __name__=="__main__":

    #load ibis wavelengths
    ibis_wavls=np.genfromtxt("ibis_wavelengths.txt")    

    #test what a typical value of SIF should be
    #fluor=0.002/convert_units_Qcm2_to_Wm2(760.78,1)
    #print(fluor)


    #test retrievals with flat spectra for
    #fluorescence and toc reflectance
    u=uvspec()
    u.options["fluorescence"]="fluorescence 0.0"
    u.options["albedo"]="albedo 1.0"
    u.options["zout"]="0.01"
    white_ref=libRadSpectra(u.run())
    white_ref.resample_to_ibis(ibis_wavls)
    white_ref.convert_units_Qcm2_to_Wm2()


    #u.options["fluorescence"]="fluorescence 0.0"
    #u.options["albedo"]="albedo 0.0"
    #black_ref=libRadSpectra(u.run())
    #black_ref.resample_to_ibis(ibis_wavls)
    #black_ref.convert_units_Qcm2_to_Wm2()


    #print(white_ref.wavl_of_min_over_band((759,770)))
    
    #flat fluorescence over a dark-ish background
    u.options["fluorescence"]="fluorescence 765440000000.0" #approx 2.0mW at 671
    u.options["albedo"]="albedo 0.1"
    flat_sif=libRadSpectra(u.run())
    flat_sif.resample_to_ibis(ibis_wavls)
    flat_sif.convert_units_Qcm2_to_Wm2()

    #example of reading a different column
    #u.run() will not run uvspec and just return
    #the name of the file containing the outputs
    flat_sif_edir=libRadSpectra(u.run(),dataCol=2)
    flat_sif_edir.resample_to_ibis(ibis_wavls)
    flat_sif_edir.convert_units_Qcm2_to_Wm2()
    
    #repeat above with a zero albedo background
    u.options["albedo"]="albedo 0.0"
    flat_sif_norefl=libRadSpectra(u.run())
    flat_sif_norefl.resample_to_ibis(ibis_wavls)
    flat_sif_norefl.convert_units_Qcm2_to_Wm2()

    #dark albedo, no SIF 
    u.options["albedo"]="albedo 0.1"
    u.options["fluorescence"]="fluorescence 0.0" 
    flat_nosif=libRadSpectra(u.run())
    flat_nosif.resample_to_ibis(ibis_wavls)
    flat_nosif.convert_units_Qcm2_to_Wm2()    

    print("outside diff:", (flat_sif.min_over_band((740,741))-flat_nosif.min_over_band((740,741)))*1000)
    print("outside diff:", (flat_sif.min_over_band((686,686.5))-flat_nosif.min_over_band((686,686.5)))*1000)
    
    #flat_sif.data*=np.pi/2.

    print("sFLD O2a: ",sif_sFLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("3FLD O2a: ",sif_3FLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("sFLD O2b: ",sif_sFLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)
    print("3FLD O2b: ",sif_3FLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)
   
    plt.plot(white_ref.wavl,white_ref.data,label="????")
    #plt.plot(flat_sif.wavl,flat_sif.data,label="SIF alb=0.1")
    #plt.plot(flat_sif_norefl.wavl,flat_sif_norefl.data,label="SIF alb=0.0")
    #plt.plot(flat_sif_edir.wavl,flat_sif_edir.data,label="direct downward flux at z")
    #plt.legend()
    #plt.show()

    
    
