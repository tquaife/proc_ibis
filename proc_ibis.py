import pickle
from copy import copy
import numpy as np

from arsf_envi_reader import envi_header as envhdr
from arsf_envi_reader.numpy_bin_reader import BilReader
from spectra_tools.spectra import Spectra
#from sif_retrieval import sif_sFLD, sif_3FLD, sif_linReg, sif_opts_FLD_Cendrero19_O2a, sif_opts_FLD_Cendrero19_O2b
from sif_retrieval import *

import sys

class ibisBilReader(BilReader):
    
    def __init__(self, ibis_filename):
        """Extends the ARSF BilReader class to add 
        elements specific to the IBIS data
        """
        super(BilReader, self).__init__(ibis_filename)

        #add wavelengths to data :
        wvl=envhdr.read_hdr_file(self.header_file)["wavelength"]
        self.wavelength=np.array(wvl.split(",")).astype(float)
        
        #scale factor to convert 
        #units from nW/cm2/sr/nm
        #to mW/m2/sr/nm
        self.units_to_mw_m2=0.01    

    
def ibis_ndvi( ibis_spect, wref, red_wband=(670, 680), nir_wband=(770, 780)):
    """Calculates a rough NDVI type index using either end 
    of the IBIS spectral range. 
    """
    red=ibis_spect.avg_over_band(red_wband)/wref.avg_over_band(red_wband)
    nir=ibis_spect.avg_over_band(nir_wband)/wref.avg_over_band(nir_wband)
    return (nir-red)/(red+nir)

if __name__=="__main__":

    #create ibis data object
    ibis_filename="./data_in/i200073b_mapped_utm30n.bil"
    ibis_data=ibisBilReader(ibis_filename)

    #read in reference data
    wref=Spectra(fname="./data_in/white_reference.csv",ftype="CSV",hdrLines=0)

    #spectra object to hold ibis data
    ibis_spect=Spectra()
    ibis_spect.wavl=copy(ibis_data.wavelength)
    
    #blank array to hold ouput data
    output_sif_o2a=np.ones([ibis_data.lines,ibis_data.samples])*-9999
    output_sif_o2b=np.ones([ibis_data.lines,ibis_data.samples])*-9999
    #output_ndvi=np.ones([ibis_data.lines,ibis_data.samples])*-9999

    #restrict the O2b waveband - the longer
    #wavelengths change too much (this may
    #also be overcome with the ridge regression)
    sif_opts_FLD_Cendrero19_O2b.in_wband=(682, 690) 

    #precompute pseudo inverse for regression
    #pseudo_o2a=sif_linReg_compute_pseudo_inverse(wref, sif_opts_FLD_Cendrero19_O2a)
    #pseudo_o2b=sif_linReg_compute_pseudo_inverse(wref, sif_opts_FLD_Cendrero19_O2b)
    #pseudo_o2a=sif_ridgeReg_compute_pseudo_inverse(wref, order=2, gamma=10000000., sif_opts=sif_opts_FLD_Cendrero19_O2a)
    pseudo_o2b=sif_ridgeReg_compute_pseudo_inverse(wref, order=2, gamma=10000000., sif_opts=sif_opts_FLD_Cendrero19_O2b)

    #loop over the IBIS data, exploiting the
    #fact that the BilReader class creates an
    #iterable object
    for (i,line) in enumerate(ibis_data):
        print(i,end=" ",flush=True)
        for j in range(ibis_data.samples):
            ibis_spect.data=copy(line[:,j])
            if not ibis_spect.data.any():
               continue
            #output_sif_o2a[i,j]=sif_linReg(ibis_spect, wref, sif_opts_FLD_Cendrero19_O2a)
            #output_sif_o2b[i,j]=sif_linReg(ibis_spect, wref, sif_opts_FLD_Cendrero19_O2b)
            #output_sif_o2a[i,j]=sif_linReg_use_precomp_inverse(ibis_spect, pseudo_o2a, sif_opts=sif_opts_FLD_Cendrero19_O2a)
            #output_sif_o2b[i,j]=sif_linReg_use_precomp_inverse(ibis_spect, pseudo_o2b, sif_opts=sif_opts_FLD_Cendrero19_O2b)
            #output_sif_o2a[i,j]=sif_ridgeReg_use_precomp_inverse(ibis_spect, pseudo_o2a, sif_opts=sif_opts_FLD_Cendrero19_O2a,out_wvl=761.0)
            output_sif_o2b[i,j]=sif_ridgeReg_use_precomp_inverse(ibis_spect, pseudo_o2b, sif_opts=sif_opts_FLD_Cendrero19_O2b,out_wvl=686.7)
            #output_ndvi[i,j]=ibis_ndvi(ibis_spect, wref)
    
    #convert units
    output_sif_o2a*=ibis_data.units_to_mw_m2
    output_sif_o2b*=ibis_data.units_to_mw_m2


    #with open('./data_out/sif_o2a_ridgeReg2.pickle', 'wb') as f:
    #    pickle.dump(output_sif_o2a, f, pickle.HIGHEST_PROTOCOL)
    with open('./data_out/sif_o2b_ridgeReg2.pickle', 'wb') as f:
        pickle.dump(output_sif_o2b, f, pickle.HIGHEST_PROTOCOL)

    #with open('./data_out/sif_o2a_linReg.pickle', 'wb') as f:
    #    pickle.dump(output_sif_o2a, f, pickle.HIGHEST_PROTOCOL)
    #with open('./data_out/sif_o2b_linReg.pickle', 'wb') as f:
    #    pickle.dump(output_sif_o2b, f, pickle.HIGHEST_PROTOCOL)
    #with open('./data_out/sif_o2a_3FLD.pickle', 'wb') as f:
    #    pickle.dump(output_sif_o2a, f, pickle.HIGHEST_PROTOCOL)
    #with open('./data_out/sif_o2b_3FLD.pickle', 'wb') as f:
    #    pickle.dump(output_sif_o2b, f, pickle.HIGHEST_PROTOCOL)
    #with open('./data_out/ndvi.pickle', 'wb') as f:
    #    pickle.dump(output_ndvi, f, pickle.HIGHEST_PROTOCOL)

    

