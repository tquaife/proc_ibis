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

    
def ibis_iterator(ibis_data,verbose=True):
    """loop across an ibisBilReader object and
    return the coordinates of the pixel (i,j)
    and a Spectra object containing the observed
    data. 
    """

    #spectra object to hold ibis data
    ibis_spect=Spectra()
    ibis_spect.wavl=copy(ibis_data.wavelength)

    #reset these counters and file position 
    #or it's only possible to call the iterator 
    #once each time a data set is opened.
    ibis_data.current_line = -1
    ibis_data.current_band = -1
    ibis_data.file_handler.seek(0)

    #loop structure
    for (i,line) in enumerate(ibis_data):
        if verbose:
            print(i,end=" ",flush=True)
        for j in range(ibis_data.samples):
            ibis_spect.data=copy(line[:,j])
            if not ibis_spect.data.any():
               continue
            yield (i,j,ibis_spect)


def ibis_sif_rFLD(wref,ibis_data,diffR,gR,diffF,gF,sif_opts,out_fn,out_wvl):
    """Apply the rFLD method to an IBIS image. Writes the results to a pickle.
    
    wref:       (Spectra object) the reference spectrum
    ibis_data:  (ibisBilReader object) the ibis data
    diffR:      (int) the order of difference constraint on reflectance
    gR:         (float) the ridge parameter for the reflectance 
    diffF:      (int) the order of difference constraint on fluorescence
    gF:         (float) the ridge parameter for the fluorescence
    sif_opts:   (sif_opts object) containing the wavelength ranges
    out_fn:     (string) the name of the file to write
    out_wvl:    (float) the wavelength of the fluorescence to output
    
    Returns:    None
    """ 
 
    #blank array to hold output data
    output_sif=np.ones([ibis_data.lines,ibis_data.samples])*-9999

    #precompute pseudo inverse for regression
    pseudo=sif_ridgeReg_compute_pseudo_inverse_RFopts(wref,\
               orderR=diffR,gammaR=gR,orderF=diffF,gammaF=gF,\
               sif_opts=sif_opts)

    #loop over image data
    for (i,j,ibis_spect) in ibis_iterator(ibis_data):
        output_sif[i,j]=sif_ridgeReg_use_precomp_inverse(ibis_spect,\
                pseudo, sif_opts=sif_opts,out_wvl=out_wvl)

    #convert units
    output_sif*=ibis_data.units_to_mw_m2

    #write data
    with open(out_fn, 'wb') as f:
        pickle.dump(output_sif, f, pickle.HIGHEST_PROTOCOL)


def ibis_sif_genericFLD(fld_func,wref,ibis_data,sif_opts,out_fn): 
    """Apply the generic method to an IBIS image. Writes the results 
    to a pickle. the FLD function must take arguments of 
    (ibis_spect, wref, sif_opts) and not require any additional
    information.
    
    fld_func:   (function) the FLD function to apply to IBIS data
    wref:       (Spectra object) the reference spectrum
    ibis_data:  (ibisBilReader object) the ibis data
    sif_opts:   (sif_opts object) containing the wavelength ranges
    out_fn:     (string) the name of the file to write
    
    Returns:    None
    """ 
 
    #blank array to hold output data
    output_sif=np.ones([ibis_data.lines,ibis_data.samples])*-9999

    #loop over image data
    for (i,j,ibis_spect) in ibis_iterator(ibis_data):
        output_sif[i,j]=fld_func(ibis_spect, wref, sif_opts=sif_opts)
        
    #convert units
    output_sif*=ibis_data.units_to_mw_m2

    #write data
    with open(out_fn, 'wb') as f:
        pickle.dump(output_sif, f, pickle.HIGHEST_PROTOCOL)


def ibis_ndvi( ibis_spect, wref, red_wband=(670, 680), nir_wband=(770, 780)):
    """Calculates a rough NDVI type index using either end 
    of the IBIS spectral range. 
    """
    red=ibis_spect.avg_over_band(red_wband)/wref.avg_over_band(red_wband)
    nir=ibis_spect.avg_over_band(nir_wband)/wref.avg_over_band(nir_wband)
    return (nir-red)/(red+nir)


if __name__=="__main__":

    #open ibis data 
    ibis_filename="./data_in/i200073b_mapped_utm30n.bil"
    ibis_data=ibisBilReader(ibis_filename)

    #read in reference data
    wref=Spectra(fname="./data_in/white_reference.csv",ftype="CSV",hdrLines=0)

    #rFLD O2a
    out_fn='./data_out/sif_o2a_rFLD.pickle'
    sif_opts=sif_opts_FLD_Cendrero19_O2a
    #ibis_sif_rFLD(wref,ibis_data,2,10E6,2,10E6,sif_opts,out_fn,761.0) 

    #rFLD O2b
    out_fn='./data_out/sif_o2b_rFLD.pickle'
    sif_opts=sif_opts_FLD_Cendrero19_O2b
    #ibis_sif_rFLD(wref,ibis_data,3,10E6,2,10E6,sif_opts,out_fn,686.7) 

    #sFLD O2a
    out_fn='./data_out/sif_o2a_sFLD.pickle'
    sif_opts=sif_opts_FLD_Cendrero19_O2a
    #ibis_sif_genericFLD(sif_sFLD_full_search,wref,ibis_data,sif_opts,out_fn) 

    #sFLD O2b
    out_fn='./data_out/sif_o2b_sFLD.pickle'
    sif_opts=sif_opts_FLD_Cendrero19_O2b
    #ibis_sif_genericFLD(sif_sFLD_full_search,wref,ibis_data,sif_opts,out_fn) 

    #3FLD O2a
    out_fn='./data_out/sif_o2a_3FLD.pickle'
    sif_opts=sif_opts_FLD_Cendrero19_O2a
    #ibis_sif_genericFLD(sif_3FLD_full_search,wref,ibis_data,sif_opts,out_fn) 

    #3FLD O2b
    out_fn='./data_out/sif_o2b_3FLD.pickle'
    sif_opts=sif_opts_FLD_Cendrero19_O2b
    #ibis_sif_genericFLD(sif_3FLD_full_search,wref,ibis_data,sif_opts,out_fn) 

    #iFLD O2a
    out_fn='./data_out/sif_o2a_iFLD.pickle'
    sif_opts=sif_opts_FLD_Cendrero19_O2a
    #ibis_sif_genericFLD(sif_iFLD_full_search,wref,ibis_data,sif_opts,out_fn) 

    #iFLD O2b
    out_fn='./data_out/sif_o2b_iFLD.pickle'
    sif_opts=sif_opts_FLD_Cendrero19_O2b
    #ibis_sif_genericFLD(sif_iFLD_full_search,wref,ibis_data,sif_opts,out_fn) 

    #The following are examples of applying the "fast" version
    #of the sFLD (a similar function exists for 3FLD). They do not 
    #currently produce great results but it demonstrates the flexibility 
    #of the ibis_sif_genericFLD() function.

    if False:
        #sFLD O2a fast algorithm    
        out_fn='./data_out/sif_o2a_sFLD_fast.pickle'
        targ=Spectra(fname="./data_in/test_spectra.csv",ftype="CSV",hdrLines=0)
        sif_opts=get_FLD_wavls_and_idx(targ, wref, sif_opts=sif_opts_FLD_Cendrero19_O2a)
        ibis_sif_genericFLD(sif_sFLD_use_predetermined_idx,wref,ibis_data,sif_opts,out_fn) 

    if False:
        #sFLD O2b fast algorithm    
        out_fn='./data_out/sif_o2b_sFLD_fast.pickle'
        targ=Spectra(fname="./data_in/test_spectra.csv",ftype="CSV",hdrLines=0)
        sif_opts=get_FLD_wavls_and_idx(targ, wref, sif_opts=sif_opts_FLD_Cendrero19_O2b)
        ibis_sif_genericFLD(sif_sFLD_use_predetermined_idx,wref,ibis_data,sif_opts,out_fn) 


