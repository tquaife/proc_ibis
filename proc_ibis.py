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


def sif_FLD( ibis_spect, wref, out_wband=(755, 757), in_wband=(761, 761.5)):
    """calculate SIF using the most basic method, taking averages 
    over the inner and outer parts of the absorption feature
    
    "in" and "out" refer to the wavelengths to be used 
    inside and outside the O2a absorption feature
    """
    Eout=wref.avg_over_band(out_wband)
    Ein=wref.avg_over_band(in_wband)
    Lout=ibis_spect.avg_over_band(out_wband)
    Lin=ibis_spect.avg_over_band(in_wband)
    return (Eout*Lin-Ein*Lout)/(Eout-Ein)


def sif_sFLD( ibis_spect, wref, band="O2a" ):
    """calculate SIF using the sFLD method defined in 
    Pilar et al. (2019) https://doi.org/10.3390/rs11080962
    
    Find the minimum and maximum values in predefined bands
    and carry out and FLD calculation with those values.
    Band definitions are from Pilar et al. (ibid).
    """
    
    if band=="O2a":
        out_wband=(745, 759) 
        in_wband=(759, 770)
    elif band=="O2b":
        out_wband=(680, 686) 
        in_wband=(686, 697)
    else:
        raise Exception("undefined absorption feature")
    
    Eout=wref.max_over_band(out_wband)
    Ein=wref.min_over_band(in_wband)
    Lout=ibis_spect.max_over_band(out_wband)
    Lin=ibis_spect.min_over_band(in_wband)
    return (Eout*Lin-Ein*Lout)/(Eout-Ein)

def sif_sFLD__test2( ibis_spect, wref, band="O2a" ):
    """calculate SIF using the sFLD method defined in 
    Pilar et al. (2019) https://doi.org/10.3390/rs11080962
    
    Find the minimum and maximum values in predefined bands
    and carry out and FLD calculation with those values.
    Band definitions are from Pilar et al. (ibid).
    """
    
    if band=="O2a":
        out_wband=(752.5, 753) 
        in_wband=(760.5, 761.0)
    elif band=="O2b":
        out_wband=(680, 686) 
        in_wband=(686, 697)
    else:
        raise Exception("undefined absorption feature")
    
    Eout=wref.avg_over_band(out_wband)
    Ein=wref.avg_over_band(in_wband)
    Lout=ibis_spect.avg_over_band(out_wband)
    Lin=ibis_spect.avg_over_band(in_wband)
    return (Eout*Lin-Ein*Lout)/(Eout-Ein)


def sif_sFLD__test( ibis_spect, wref, band="O2a" ):
    """calculate SIF using the sFLD method defined in 
    Pilar et al. (2019) https://doi.org/10.3390/rs11080962
    
    Find the minimum and maximum values in predefined bands
    and carry out and FLD calculation with those values.
    Band definitions are from Pilar et al. (ibid).
    """
    
    if band=="O2a":
        out_wl_left=753  
        #out_wl_right=771 
        in_wl=760
    elif band=="O2b":
        out_wl_left=686.5 
        #out_wl_right=697.5 
        in_wl=687.0
    else:
        raise Exception("undefined absorption feature")
    
    (_,Eout)=wref.closest_to_wavl(out_wl_left)
    (_,Ein)=wref.closest_to_wavl(in_wl)
    (_,Lout)=ibis_spect.closest_to_wavl(out_wl_left)
    (_,Lin)=ibis_spect.closest_to_wavl(in_wl)
    return (Eout*Lin-Ein*Lout)/(Eout-Ein)





def sif_3FLD( ibis_spect, wref, band="O2a" ):
    """calculate SIF using the 3FLD method defined in 
    Damm et al. (2011) https://doi.org/10.1016/j.rse.2011.03.011
    with band definitions from
    Pilar et al. (2019) https://doi.org/10.3390/rs11080962
    
    Find the minimum and maximum values in predefined bands
    and carry out and 3FLD calculation with those values.
    """
    
    if band=="O2a":
        out_wband_left=(745, 759) 
        out_wband_right=(770, 780) 
        in_wband=(759, 770)
    elif band=="O2b":
        out_wband_left=(680, 686) 
        out_wband_right=(697, 698) 
        in_wband=(686, 697)
    else:
        raise Exception("undefined absorption feature")
    
    Eout_L=wref.max_over_band(out_wband_left)
    Eout_R=wref.max_over_band(out_wband_right)
    Ein=wref.min_over_band(in_wband)
    
    Lout_L=ibis_spect.max_over_band(out_wband_left)
    Lout_R=ibis_spect.max_over_band(out_wband_right)
    Lin=ibis_spect.min_over_band(in_wband)

    Eout_L_wvl=wref.wavl_of_max_over_band(out_wband_left)
    Eout_R_wvl=wref.wavl_of_max_over_band(out_wband_right)

    Lin_wvl=ibis_spect.wavl_of_min_over_band(in_wband)

    w21=(Eout_R_wvl-Lin_wvl)/(Eout_R_wvl-Eout_L_wvl)
    w22=(Lin_wvl-Eout_L_wvl)/(Eout_R_wvl-Eout_L_wvl)

    Ecorr=Ein/(w21*Eout_L+w22*Eout_R)

    return (Lin-Ecorr*(w21*Lout_L+w22*Lout_R))/(1-Ecorr)


def sif_3FLD__test( ibis_spect, wref, band="O2a" ):
    """
    TEST VERSION - DO NOT USE
    
    calculate SIF using the 3FLD method defined in 
    Damm et al. (2011) https://doi.org/10.1016/j.rse.2011.03.011
    with band definitions from
    Pilar et al. (2019) https://doi.org/10.3390/rs11080962
    
    Find the minimum and maximum values in predefined bands
    and carry out and 3FLD calculation with those values.
    """
    
    if band=="O2a":
        out_wl_left=753  
        out_wl_right=771 
        in_wl=760
    elif band=="O2b":
        out_wl_left=686.5 
        out_wl_right=697.5 
        in_wl=687.0
    else:
        raise Exception("undefined absorption feature")
    
    (Eout_L_wvl,Eout_L)=wref.closest_to_wavl(out_wl_left)          
    (Eout_R_wvl,Eout_R)=wref.closest_to_wavl(out_wl_right)          
    (_,Ein)=wref.closest_to_wavl(in_wl)

    (Lout_L_wvl,Lout_L)=ibis_spect.closest_to_wavl(out_wl_left)          
    (Lout_R_wvl,Lout_R)=ibis_spect.closest_to_wavl(out_wl_right)          
    (Lin_wvl,Lin)=ibis_spect.closest_to_wavl(in_wl)

    w21=(Eout_R_wvl-Lin_wvl)/(Eout_R_wvl-Eout_L_wvl)
    w22=(Lin_wvl-Eout_L_wvl)/(Eout_R_wvl-Eout_L_wvl)

    Ecorr=Ein/(w21*Eout_L+w22*Eout_R)

    return (Lin-Ecorr*(w21*Lout_L+w22*Lout_R))/(1-Ecorr)




    
def ibis_ndvi( ibis_spect, wref, red_wband=(670, 680), nir_wband=(770, 780)):
    """Calculates a rough NDVI type index using either end 
    of the IBIS spectral range. 
    """
    red=ibis_spect.avg_over_band(red_wband)/wref.avg_over_band(red_wband)
    nir=ibis_spect.avg_over_band(nir_wband)/wref.avg_over_band(nir_wband)
    return (nir-red)/(red+nir)


if __name__=="__main__":

    #create ibis data object
    ibis_filename="./data/i200073b_mapped_utm30n.bil"
    ibis_data=ibisBilReader(ibis_filename)

    #read in reference data
    wref=Spectra(fname="data/white_reference.csv",ftype="CSV",hdrLines=0)

    #spectra object to hold ibis data
    ibis_spect=Spectra()
    ibis_spect.wavl=copy(ibis_data.wavelength)
    
    #blank array to hold ouput data
    output_sif_o2a=np.ones([ibis_data.lines,ibis_data.samples])*-9999
    output_sif_o2b=np.ones([ibis_data.lines,ibis_data.samples])*-9999
    #output_ndvi=np.ones([ibis_data.lines,ibis_data.samples])*-9999

    #loop over the IBIS data, exploiting the
    #fact that the BilReader class creates an
    #iterable object
    for (i,line) in enumerate(ibis_data):
        for j in range(ibis_data.samples):
            ibis_spect.data=copy(line[:,j])
            if not ibis_spect.data.any():
               continue
            output_sif_o2a[i,j]=sif_3FLD__test(ibis_spect, wref, band="O2a")
            output_sif_o2b[i,j]=sif_3FLD__test(ibis_spect, wref, band="O2b")
            #output_ndvi[i,j]=ibis_ndvi(ibis_spect, wref)
    
    #convert units
    output_sif_o2a*=ibis_data.units_to_mw_m2
    output_sif_o2b*=ibis_data.units_to_mw_m2

    with open('sif_o2a_3FLD__test.pickle', 'wb') as f:
        pickle.dump(output_sif_o2a, f, pickle.HIGHEST_PROTOCOL)
    with open('sif_o2b_3FLD__test.pickle', 'wb') as f:
        pickle.dump(output_sif_o2b, f, pickle.HIGHEST_PROTOCOL)
    #with open('ndvi.pickle', 'wb') as f:
    #    pickle.dump(output_ndvi, f, pickle.HIGHEST_PROTOCOL)

    

