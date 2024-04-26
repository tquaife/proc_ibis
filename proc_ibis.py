import pickle
from copy import copy
import numpy as np

from arsf_envi_reader import envi_header as envhdr
from arsf_envi_reader.numpy_bin_reader import BilReader
from spectra_tools.spectra import Spectra

import sys

class ibisBilReader(BilReader):
    
    def __init__(self, ibis_filename ):
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


def sif_FLD( ibis_spect, reference, out_wband=(755, 757), in_wband=(761, 761.5)):
    """calculate SIF using the most basic method
    
    "in" and "out" refer to the wavelengths to be used 
    inside and outside the O2a absorption feature
    """
    Eout=wref.avg_over_band(out_wband)
    Ein=wref.avg_over_band(in_wband)
    Lout=ibis_spect.avg_over_band(out_wband)
    Lin=ibis_spect.avg_over_band(in_wband)
    return (Eout*Lin-Ein*Lout)/(Eout-Ein)

    

if __name__=="__main__":

    #create ibis data object
    ibis_filename="./data/i200073b_mapped_utm30n.bil"
    ibis_data=ibisBilReader(ibis_filename)

    #read in reference data
    wref=Spectra(fname="white_reference.csv",ftype="CSV",hdrLines=0)

    #spectra object to hold ibis data
    ibis_spect=Spectra()
    ibis_spect.wavl=copy(ibis_data.wavelength)
    
    #blank array to hold SIF data
    output_sif=np.ones([ibis_data.lines,ibis_data.samples])*-9999

    #loop over the IBIS data, exploiting the
    #fact that the BilReader class creates an
    #iterable object
    for (i,line) in enumerate(ibis_data):
        for j in range(ibis_data.samples):
            ibis_spect.data=copy(line[:,j])
            if not ibis_spect.data.any():
               continue
            output_sif[i,j]=sif_FLD(ibis_spect, wref)
    
    #convert units
    output_sif*=ibis_data.units_to_mw_m2

    with open('sif.pickle', 'wb') as f:
        pickle.dump(output_sif, f, pickle.HIGHEST_PROTOCOL)



