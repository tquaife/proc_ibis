import sys
sys.path.append("../")

import numpy as np
from matplotlib import pyplot as plt

from ibis_simulator import libRadSpectra, convert_units_Qcm2_to_Wm2
from sif_retrieval import sif_sFLD, sif_3FLD, sif_iFLD, sif_opts_FLD_Cendrero19_O2a, sif_opts_FLD_Cendrero19_O2b
from run_libradtran import uvspec


if __name__=="__main__":

    #setting this to False means libRadTran will
    #only run if the output file for the run doesn't exist
    overwrite=True
    #overwrite=False

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
    u.options["zout"]="1.0"
    #u.options["zout"]="sur"
    white_ref=libRadSpectra(u.run(overwrite=overwrite))
    white_ref.convert_units_Qcm2_to_Wm2()
    white_ref.resample_to_ibis(ibis_wavls)

    white_ref_edir=libRadSpectra(u.run(overwrite=overwrite),dataCol=2)
    white_ref_edir.convert_units_Qcm2_to_Wm2()
    white_ref_edir.resample_to_ibis(ibis_wavls)

    white_ref_edif=libRadSpectra(u.run(overwrite=overwrite),dataCol=3)
    white_ref_edif.convert_units_Qcm2_to_Wm2()
    white_ref_edif.resample_to_ibis(ibis_wavls)

    white_ref_eup=libRadSpectra(u.run(overwrite=overwrite),dataCol=4)
    white_ref_eup.convert_units_Qcm2_to_Wm2()
    white_ref_eup.resample_to_ibis(ibis_wavls)
    
    #flat fluorescence over a dark-ish background
    u.options["fluorescence"]="fluorescence 765440000000.0" #approx 2.0mW at 671
    u.options["albedo"]="albedo 0.1"
    flat_sif=libRadSpectra(u.run(overwrite=overwrite))
    flat_sif.resample_to_ibis(ibis_wavls)
    flat_sif.convert_units_Qcm2_to_Wm2()

    flat_sif_edir=libRadSpectra(u.run(overwrite=overwrite),dataCol=2)
    flat_sif_edir.convert_units_Qcm2_to_Wm2()
    flat_sif_edir.resample_to_ibis(ibis_wavls)

    flat_sif_edif=libRadSpectra(u.run(overwrite=overwrite),dataCol=3)
    flat_sif_edif.convert_units_Qcm2_to_Wm2()
    flat_sif_edif.resample_to_ibis(ibis_wavls)
    
    #repeat above with a zero albedo background
    u.options["albedo"]="albedo 0.0"
    flat_sif_norefl=libRadSpectra(u.run(overwrite=overwrite))
    flat_sif_norefl.resample_to_ibis(ibis_wavls)
    flat_sif_norefl.convert_units_Qcm2_to_Wm2()

    #dark albedo, no SIF 
    u.options["albedo"]="albedo 0.1"
    u.options["fluorescence"]="fluorescence 0.0" 
    flat_nosif=libRadSpectra(u.run(overwrite=overwrite))
    flat_nosif.resample_to_ibis(ibis_wavls)
    flat_nosif.convert_units_Qcm2_to_Wm2()    


    flat_nosif_edir=libRadSpectra(u.run(overwrite=overwrite),dataCol=2)
    flat_nosif_edir.convert_units_Qcm2_to_Wm2()
    flat_nosif_edir.resample_to_ibis(ibis_wavls)

    flat_nosif_edif=libRadSpectra(u.run(overwrite=overwrite),dataCol=3)
    flat_nosif_edif.convert_units_Qcm2_to_Wm2()
    flat_nosif_edif.resample_to_ibis(ibis_wavls)


    #print("outside diff:", (flat_sif.avg_over_band((740,741))-flat_nosif.avg_over_band((740,741)))*1000)
    #print("outside diff:", (flat_sif.avg_over_band((686,686.5))-flat_nosif.avg_over_band((686,686.5)))*1000)
    print("diff @ 761:", (flat_sif.avg_over_band((761,761.5))-flat_nosif.avg_over_band((761,761.5)))*1000)
    
    #flat_sif.data*=np.pi/2.

    print("sFLD O2a: ",sif_sFLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("sFLD O2b: ",sif_sFLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)
    print("3FLD O2a: ",sif_3FLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("3FLD O2b: ",sif_3FLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)
    print("iFLD O2a: ",sif_iFLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("iFLD O2b: ",sif_iFLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)
   
    
    plt.plot(flat_sif.wavl,flat_sif.data,label="SIF alb=0.1")
    plt.plot(flat_sif_norefl.wavl,flat_sif_norefl.data,label="SIF alb=0.0")
    
    #these two are the same at the surface
    plt.plot(flat_sif_edir.wavl,flat_sif_edir.data+flat_sif_edif.data,"k-",label="total downward flux at z (flat sif)")
    plt.plot(flat_sif_edir.wavl,flat_nosif_edir.data+flat_nosif_edif.data,"--",label="total downward flux at z (flat NO sif)")

    #---
    #presumably the difference between the ones above and below this comment
    #is multiple scattering between the surface an the atmosphere
    #---

    #these are all the same at the surface:
    plt.plot(white_ref.wavl,white_ref.data*np.pi,"b-",label="white ref * pi",linewidth=3.)    
    plt.plot(white_ref_edir.wavl,white_ref_edir.data+white_ref_edif.data,"r--",label="total downward flux at z (white ref)",linewidth=2.)
    plt.plot(white_ref_eup.wavl,white_ref_eup.data,"y-.",label="eup from white_ref")

    plt.title("Z out ="+u.options["zout"]+"km")

    plt.legend()
    plt.show()

    
    
