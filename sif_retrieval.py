import numpy as np
from spectra_tools.spectra import Spectra

class sif_opts_base: pass

# Band definitions from Pilar et al. (2019) 
# https://doi.org/10.3390/rs11080962
sif_opts_FLD_Cendrero19_O2a=sif_opts_base()
sif_opts_FLD_Cendrero19_O2a.in_wband=(759, 770) 
sif_opts_FLD_Cendrero19_O2a.outL_wband=(745, 759)
sif_opts_FLD_Cendrero19_O2a.outR_wband=(770, 780) 
sif_opts_FLD_Cendrero19_O2b=sif_opts_base()
sif_opts_FLD_Cendrero19_O2b.in_wband=(686, 697) 
sif_opts_FLD_Cendrero19_O2b.outL_wband=(680, 686)  
sif_opts_FLD_Cendrero19_O2b.outR_wband=(697, 698) 


def sif_sFLD( target_spect, wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a ):
    """calculate SIF using the sFLD method defined in 
    Cendrero-Mateo et al. (2019) https://doi.org/10.3390/rs11080962
    """
    #get the L and E at the local maxima (in E) that
    #is closest to the absorption feature    
    tmp_spect=get_local_maxima(wref_spect,sif_opts.outL_wband[0],sif_opts.outL_wband[1])
    Eout=tmp_spect.data[-1]    
    tmp_spect=get_local_maxima(target_spect,sif_opts.outL_wband[0],sif_opts.outL_wband[1])
    Lout=tmp_spect.data[-1]
    #use this to get Lout at the same wavl as Eout:    
    #(_,Lout)=target_spect.closest_to_wavl(tmp_spect.wavl[-1])
    
    #get L and E in the absorption feature
    #assumes position of minima is same in L and E
    Ein=wref_spect.min_over_band(sif_opts.in_wband)
    Lin=target_spect.min_over_band(sif_opts.in_wband)

    return (Eout*Lin-Ein*Lout)/(Eout-Ein)


def sif_3FLD( target_spect, wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a  ):
    """calculate SIF using the 3FLD method defined in 
    Damm et al. (2011) https://doi.org/10.1016/j.rse.2011.03.011
    """

    Ein=wref_spect.min_over_band(sif_opts.in_wband)
    tmp_spect=get_local_maxima(wref_spect,sif_opts.outL_wband[0],sif_opts.outL_wband[1])    
    Eout_L=tmp_spect.data[-1]
    Eout_L_wvl=tmp_spect.wavl[-1]
    tmp_spect=get_local_maxima(wref_spect,sif_opts.outR_wband[0],sif_opts.outR_wband[1])    
    Eout_R=tmp_spect.data[0]
    Eout_R_wvl=tmp_spect.wavl[0]
    
    Lin=target_spect.min_over_band(sif_opts.in_wband)
    Lin_wvl=target_spect.wavl_of_min_over_band(sif_opts.in_wband)
    tmp_spect=get_local_maxima(target_spect,sif_opts.outL_wband[0],sif_opts.outL_wband[1])    
    Lout_L=tmp_spect.data[-1]
    tmp_spect=get_local_maxima(target_spect,sif_opts.outR_wband[0],sif_opts.outR_wband[1])    
    Lout_R=tmp_spect.data[0]    

    w21=(Eout_R_wvl-Lin_wvl)/(Eout_R_wvl-Eout_L_wvl)
    w22=(Lin_wvl-Eout_L_wvl)/(Eout_R_wvl-Eout_L_wvl)
    Ecorr=Ein/(w21*Eout_L+w22*Eout_R)

    return (Lin-Ecorr*(w21*Lout_L+w22*Lout_R))/(1-Ecorr)


def get_local_minima(in_spect, wl_beg, wl_end):
    """return a Spectra object that only contains
    data where the there is a local minima between
    wl_beg and wl_end.
    
    Works by inverting the data in in_spect and
    then calling get_local_maxima(). Needs to
    then un-invert the data. (it's a bit messy to do
    it this way, but if we need to change the 
    get_local_maxima() function we don't need to 
    worry about translating the changes to this one.)
    """
    in_spect.data*=-1

    out_spect=get_local_maxima(in_spect,wl_beg,wl_end)
    out_spect.data*=-1
    
    in_spect.data*=-1
    return out_spect

def get_local_maxima(in_spect, wl_beg, wl_end):
    """return a Spectra object that only contains
    data where the there is a local maxima between
    wl_beg and wl_end.
    """
    idx_beg = int(np.abs(in_spect.wavl-wl_beg).argmin())
    idx_end = int(np.abs(in_spect.wavl-wl_end).argmin())

    out_spect=Spectra()
    for i in range(idx_beg,idx_end):
        if in_spect.data[i]>in_spect.data[i-1] and\
           in_spect.data[i]>in_spect.data[i+1]:
            out_spect.data=np.append(out_spect.data,in_spect.data[i])
            out_spect.wavl=np.append(out_spect.wavl,in_spect.wavl[i])

    #catch the case where there are no local maxima
    #and return the largest value
    if len(out_spect.data)==0:
        out_spect.data=np.append(out_spect.data,in_spect.max_over_band((wl_beg, wl_end)))
        out_spect.wavl=np.append(out_spect.data,in_spect.wavl_of_max_over_band((wl_beg, wl_end)))

    return out_spect




