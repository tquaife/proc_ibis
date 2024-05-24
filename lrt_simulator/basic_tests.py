import sys
sys.path.append("../")

from copy import copy, deepcopy

import numpy as np
from matplotlib import pyplot as plt

from ibis_simulator import libRadSpectra, convert_units_Qcm2_to_Wm2
from sif_retrieval import sif_sFLD, sif_3FLD, sif_opts_FLD_Cendrero19_O2a, sif_opts_FLD_Cendrero19_O2b
from run_libradtran import uvspec

def simplest_case():
    """Runs tests at surface level against simulated white 
    reference radiance with the target radiance generated in 
    Python (i.e. rather than using libratran). The target
    reflectance and fluorescence are spectrally flat. 
    
    All FLD retrieval methods should be able to get the
    values of the SIF, so this function can be used as
    a sanity check.
    """
    
    #load ibis wavelengths
    ibis_wavls=np.genfromtxt("ibis_wavelengths.txt")    

    dataCol=2
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
    print("F: %0.1f alb: %0.2f"%(F*1000,alb))
    print("sFLD O2a:",sif_sFLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("sFLD O2b:",sif_sFLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)
    print("3FLD O2a:",sif_3FLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("3FLD O2b:",sif_3FLD(flat_sif,white_ref,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)


def simplest_w_lrt_upwelling():
    
    #load ibis wavelengths
    ibis_wavls=np.genfromtxt("ibis_wavelengths.txt")    

    dataCol=1
    overwriteLrtOuts=True

    F=0.015
    #F=F/convert_units_Qcm2_to_Wm2(761,1)
    #print(convert_units_Qcm2_to_Wm2(761,F))

    alb=0.1
    
    #generate reference spectra:
    u=uvspec()
    #u.options["fluorescence"]="fluorescence %f"%0.0
    u.options["fluorescence"]="fluorescence_file fluor_flat_0p015_hires.flu"
    u.options["albedo"]="albedo %f"%alb
    u.options["zout"]="sur"
    #u.options["zout"]=0.01
    u.options["umu"]=1.0
    
    #get irradiance by adding edir and edn
    #note: all lrt runs after the first one should
    #be set to overwrite=False as there is no need
    #to re-run lrt, just need the output file name
    #(which is returned by u.run() - this could 
    #be clearer!)
    irradiance=libRadSpectra(u.run(overwrite=overwriteLrtOuts),dataCol=2)             
    tmp=libRadSpectra(u.run(overwrite=False),dataCol=3)
    irradiance.data+=tmp.data
    irradiance.resample_to_ibis(ibis_wavls)
    irradiance.convert_units_Qcm2_to_Wm2()
    irradiance.data/=np.pi

    #get the upwelling radiance:
    flat_sif=libRadSpectra(u.run(overwrite=False),dataCol=1)
    flat_sif.resample_to_ibis(ibis_wavls)
    flat_sif.convert_units_Qcm2_to_Wm2()
    #flat_sif.data+=F



    #check the inversion algorithms
    print("F: %0.1f alb: %0.2f"%(F*1000,alb))
    print("sFLD O2a:",sif_sFLD(flat_sif,irradiance,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("sFLD O2b:",sif_sFLD(flat_sif,irradiance,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)
    print("3FLD O2a:",sif_3FLD(flat_sif,irradiance,sif_opts=sif_opts_FLD_Cendrero19_O2a)*1000)
    print("3FLD O2b:",sif_3FLD(flat_sif,irradiance,sif_opts=sif_opts_FLD_Cendrero19_O2b)*1000)

    #plt.plot(irradiance.wavl,(flat_sif.data)/irradiance.data)
    #plt.plot(irradiance.wavl,irradiance.data)
    #plt.plot(flat_sif.wavl,flat_sif.data)
    #plt.show()


if __name__=="__main__":

  
    #simplest_case()
    simplest_w_lrt_upwelling()
    
    
