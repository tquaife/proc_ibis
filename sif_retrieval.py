import numpy as np

import scipy.linalg as linalg
import scipy.interpolate as interpolate

from copy import copy, deepcopy
from spectra_tools.spectra import Spectra

class sif_opts_base: pass

# Band definitions from Cendrero et al. (2019) 
# https://doi.org/10.3390/rs11080962
sif_opts_FLD_Cendrero19_O2a=sif_opts_base()
sif_opts_FLD_Cendrero19_O2a.in_wband=(759, 770) 
sif_opts_FLD_Cendrero19_O2a.outL_wband=(745, 759)
sif_opts_FLD_Cendrero19_O2a.outR_wband=(770, 780) 
sif_opts_FLD_Cendrero19_O2b=sif_opts_base()
sif_opts_FLD_Cendrero19_O2b.in_wband=(686, 697) 
sif_opts_FLD_Cendrero19_O2b.outL_wband=(680, 686)  
sif_opts_FLD_Cendrero19_O2b.outR_wband=(697, 698) 

sif_opts_rFLD_custom_O2a=sif_opts_base()
sif_opts_rFLD_custom_O2b=sif_opts_base()
sif_opts_rFLD_custom_O2a.in_wband=(759.5, 765) 
sif_opts_rFLD_custom_O2b.in_wband=(686, 697) 


def diff_matrix(size,order=1):
    """returns a difference operator of the type 
    described by Twomey in "introduction to the 
    mathematics..." (see pages 124-127)
    
    size: (int) size of matrix to be returned
    order: (int) the order of the difference matrix 
    """
    D=np.eye(size+order)
    D=np.diff(D,n=order)
    D[0:order,:]=0
    D=D[0:size,:]
    return D


def sif_ridgeReg_compute_pseudo_inverse(wref_spect, order=1, gamma=1., sif_opts=sif_opts_FLD_Cendrero19_O2a):
    """Wrapper function for sif_ridgeReg_compute_pseudo_inverse_RFopts 
    to preserve previous functionality. Passes a single value of order and 
    gamma instead of separate ones for the reflectance and fluorescence.  
    
    Very likely to be deleted in the near future
    """
    return sif_ridgeReg_compute_pseudo_inverse_RFopts(wref_spect, orderR=order, \
                          gammaR=gamma, orderF=order, gammaF=gamma,sif_opts=sif_opts)


def sif_ridgeReg_compute_pseudo_inverse_RFopts(wref_spect, orderR=1, gammaR=1., orderF=1,
                             gammaF=1.,sif_opts=sif_opts_FLD_Cendrero19_O2a,return_k=False):
    """ Calculate the pseudo inverse for the Ridge Regression
    problem.
    
    Returns the matrix: (K^TK+gB^TB)^{-1}K^T as a 2D numpy array
    
    wref_spect:   (Spectra object) the white reference spectra
    orderR:       (int) the order of difference constraint for the reflectance 
    gammaR:       (float) the magnitude of the reflectance constraint 
    orderF:       (int) the order of difference constraint for the fluorescence 
    gammaF:       (float) the magnitude of the fluorescence constraint 
    sif_opts:     (sif_opts object) containing the wavelength ranges 
    """
    wl_min=sif_opts.in_wband[0]
    wl_max=sif_opts.in_wband[1]

    wref=copy(wref_spect)    
    wref.trim(wl_min,wl_max)
    
    size=len(wref.data)
    
    #make the design matrix 
    K=np.zeros((size,size*2))
    for i in range(size):
        K[i,i]=wref.data[i]
        K[i,i+size]=1.0
    if return_k:
        return K

    KTK=np.matmul(K.T,K)
    
    #make the block matrix for the
    #ridge regression constraints
    DR=diff_matrix(size,order=orderR)
    DF=diff_matrix(size,order=orderF)
    B=linalg.block_diag(gammaR*DR,gammaF*DF)
    BTB=np.matmul(B.T,B)
    
    return(np.matmul(np.linalg.inv(KTK+BTB),K.T))


def sif_ridgeReg_use_precomp_inverse( target_spect, pseudo_inv,sif_opts=sif_opts_FLD_Cendrero19_O2a, out_wvl=None ):

    if out_wvl is None:
        out_wvl=targ.wavl[int(size/2.)]

    (R,F)=sif_ridgeReg_use_precomp_inverse_full_output(target_spect,pseudo_inv,sif_opts=sif_opts)
    return(F.closest_to_wavl(out_wvl)[1])


def sif_ridgeReg_use_precomp_inverse_full_output( target_spect, pseudo_inv,sif_opts=sif_opts_FLD_Cendrero19_O2a ):
    """ Compute the fluorescence using ridge regression. Similar in concept 
    to iFLD in that it account for change in R and F over the window,
    but here they are both determined in a single step. Allows for all
    data points from an absorption feature to be used in the retrieval. 
    
    Returns the fluorescence in the units of the input data 
    
    target_spect:    (Spectra object) observed radiance from a fluorescing target 
    pseudo_inv:      (2D numpy array) the precomputed inverse matrix
    sif_opts:        (sif_opts object) containing the wavelength ranges
    out_wvl:         (float, optional) the wavelength to compute the fluorescence at
                     if None the return F in the middle of the window  
    """
    wl_min=sif_opts.in_wband[0]
    wl_max=sif_opts.in_wband[1]

    targ=copy(target_spect)
    targ.trim(wl_min,wl_max)
    size=len(targ.data)

    F=deepcopy(targ)
    R=deepcopy(targ)

    #if out_wvl is None:
    #    out_wvl=targ.wavl[int(size/2.)]

    out=np.matmul(pseudo_inv,targ.data)
    #copy SIF into F
    F.data=out[size:]
    #copy reflectance into F
    R.data=out[:size]

    return((R,F))
    #return(targ.closest_to_wavl(out_wvl)[1])
    
    
def sif_linReg_compute_pseudo_inverse( wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a ):
    """use linear regression to solve the equation: 
    I_up_target=I_up_reference*(R_target/R_reference)+SIF
    
    This is the same as sif_linReg but uses a precomputed 
    inverse matrix, e.g. from sif_linReg_compute_pseudo_inverse()
    
    Returns the matrix: (K^TK)^{-1}K^T as a 2D numpy array
    
    wref_spect:      (Spectra object) the white reference radiance spectra
    sif_opts:        (sif_opts object) containing the wavelength ranges
    """
    wl_min=sif_opts.in_wband[0]
    wl_max=sif_opts.in_wband[1]

    wref=copy(wref_spect)    
    wref.trim(wl_min,wl_max)

    X=np.ones((len(wref.data),2))
    X[:,0]=np.array(wref.data)
    
    return(np.linalg.pinv(X))


def sif_linReg_use_precomp_inverse( target_spect, pseudo_inv, sif_opts=sif_opts_FLD_Cendrero19_O2a ):
    """use linear regression to solve the equation: 
    I_up_target=I_up_reference*(R_target/R_reference)+SIF
    
    This is the same as sif_linReg but uses a precomputed 
    inverse matrix, e.g. from sif_linReg_compute_pseudo_inverse()
    
    returns the fluorescence (assumed uniform over the window)
    
    target_spect:    (Spectra object) observed radiance from a fluorescing target 
    pseudo_inv:      (2D numpy array) the precomputed inverse matrix
    sif_opts:        (sif_opts object) containing the wavelength ranges
    """
    wl_min=sif_opts.in_wband[0]
    wl_max=sif_opts.in_wband[1]

    targ=copy(target_spect)
    targ.trim(wl_min,wl_max)

    out=np.matmul(pseudo_inv,targ.data)

    return(out[1])


def sif_linReg( target_spect, wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a ):
    """use linear regression to solve the equation: 
    I_up_target=I_up_reference*(R_target/R_reference)+SIF
    
    This is similar in concept to sFLD but uses all the data points with the 
    window. Can be quite slow if the matrix used is big but, in principle the
    inverse could be pre-computed for circumstances where the reference spectra
    is static (for example a flight campaign).
    """
    wl_min=sif_opts.in_wband[0]
    wl_max=sif_opts.in_wband[1]

    wref=copy(wref_spect)    
    wref.trim(wl_min,wl_max)
    targ=copy(target_spect)
    targ.trim(wl_min,wl_max)

    y=np.array(targ.data)
    X=np.ones((len(wref.data),2))
    X[:,0]=np.array(wref.data)

    [ans, resd, rnk, sng]=np.linalg.lstsq(X, y, rcond=None)

    return(ans[1])


def get_FLD_wavls_and_idx( target_spect, wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a ):
    """Get the wavelengths for FLD calculations. This allows them to be 
    precomputed which, for some use cases, will save a lot of processing 
    time. May not be suitable for all applications though - if in doubt 
    use the "full search" FLD functions or the ridge regression.

    target_spect:    (Spectra object) observed radiance from a fluorescing target 
    wref_spect:      (Spectra object) the white reference radiance spectra
    sif_opts:        (sif_opts object) containing the wavelength ranges

    returns:         (sif_opts object) modified to include the wavelengths
                     and array indices of relevant points on the spectrum
                     
    note: the sif_opts object is actually modified in place, so it is not
    necessary to return it. However, imho, it makes for more readable code 
    in the Python function that calls this.
    """

    #Ein
    sif_opts.Ein_wvl=wref_spect.wavl_of_min_over_band(sif_opts.in_wband)
    sif_opts.Ein_idx=wref_spect.closest_to_wavl_idx(sif_opts.Ein_wvl)
    
    #Eout_L
    tmp_spect=get_local_maxima(wref_spect,sif_opts.outL_wband[0],sif_opts.outL_wband[1])    
    sif_opts.Eout_L_wvl=tmp_spect.wavl[-1]
    sif_opts.Eout_L_idx=wref_spect.closest_to_wavl_idx(sif_opts.Eout_L_wvl)
    
    #Eout_R
    tmp_spect=get_local_maxima(wref_spect,sif_opts.outR_wband[0],sif_opts.outR_wband[1])    
    sif_opts.Eout_R_wvl=tmp_spect.wavl[0]
    sif_opts.Eout_R_idx=wref_spect.closest_to_wavl_idx(sif_opts.Eout_R_wvl)
    
    #Lin
    sif_opts.Lin_wvl=target_spect.wavl_of_min_over_band(sif_opts.in_wband)
    sif_opts.Lin_idx=target_spect.closest_to_wavl_idx(sif_opts.Lin_wvl)

    #Lout_L
    tmp_spect=get_local_maxima(target_spect,sif_opts.outL_wband[0],sif_opts.outL_wband[1])    
    sif_opts.Lout_L_wvl=tmp_spect.data[-1]
    sif_opts.Lout_L_idx=wref_spect.closest_to_wavl_idx(sif_opts.Lout_L_wvl)

    #Lout_R
    tmp_spect=get_local_maxima(target_spect,sif_opts.outR_wband[0],sif_opts.outR_wband[1])    
    sif_opts.Lout_R_wvl=tmp_spect.data[0]    
    sif_opts.Lout_R_idx=wref_spect.closest_to_wavl_idx(sif_opts.Lout_R_wvl)

    return sif_opts


def sif_sFLD_use_predetermined_idx( target_spect, wref_spect, sif_opts ):
    """calculate SIF using the sFLD method defined in 
    Cendrero-Mateo et al. (2019) https://doi.org/10.3390/rs11080962

    The sif_opts object must contain the indices of the data needed 
    to calculate the FLD. e.g. by calling get_FLD_wavelengths()
    """
    Eout=wref_spect.data[sif_opts.Eout_L_idx]    
    Lout=target_spect.data[sif_opts.Lout_L_idx]    
    Ein=wref_spect.data[sif_opts.Ein_idx]
    Lin=target_spect.data[sif_opts.Lin_idx] 

    return (Eout*Lin-Ein*Lout)/(Eout-Ein)


def sif_3FLD_use_predetermined_idx( target_spect, wref_spect, sif_opts ):
    """calculate SIF using the 3FLD method defined in 
    Damm et al. (2011) https://doi.org/10.1016/j.rse.2011.03.011

    The sif_opt objects must contain the indices of the data needed 
    to calculate the FLD. e.g. by calling get_FLD_wavelengths()
    """
    Eout_L=wref_spect.data[sif_opts.Eout_L_idx]    
    Lout_L=target_spect.data[sif_opts.Lout_L_idx]    
    Eout_R=wref_spect.data[sif_opts.Eout_R_idx]    
    Lout_R=target_spect.data[sif_opts.Lout_R_idx]    
    Ein=wref_spect.data[sif_opts.Ein_idx]
    Lin=target_spect.data[sif_opts.Lin_idx] 

    w21=(sif_opts.Eout_R_wvl-sif_opts.Lin_wvl)/(sif_opts.Eout_R_wvl-sif_opts.Eout_L_wvl)
    w22=(sif_opts.Lin_wvl-sif_opts.Eout_L_wvl)/(sif_opts.Eout_R_wvl-sif_opts.Eout_L_wvl)
    Ecorr=Ein/(w21*Eout_L+w22*Eout_R)

    return (Lin-Ecorr*(w21*Lout_L+w22*Lout_R))/(1-Ecorr)


def sif_sFLD( target_spect, wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a ):
    return sif_sFLD_full_search( target_spect, wref_spect, sif_opts=sif_opts )

def sif_sFLD_full_search( target_spect, wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a ):
    """calculate SIF using the sFLD method defined in 
    Cendrero-Mateo et al. (2019) https://doi.org/10.3390/rs11080962

    The "full_search" version finds the wavelengths one which to 
    perform the calculations itself. If processing large amounts 
    of data this can slow down the calculation significantly.
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

def sif_3FLD( target_spect, wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a ):
    return sif_3FLD_full_search( target_spect, wref_spect, sif_opts=sif_opts )

def sif_3FLD_full_search( target_spect, wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a ):
    """calculate SIF using the 3FLD method defined in 
    Damm et al. (2011) https://doi.org/10.1016/j.rse.2011.03.011
    
    The "full_search" version finds the wavelengths one which to 
    perform the calculations itself. If processing large amounts 
    of data this can slow down the calculation significantly.
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


def sif_iFLD(target_spect, wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a):
    return sif_iFLD_full_search(target_spect, wref_spect, sif_opts=sif_opts)

def sif_iFLD_full_search(target_spect, wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a):
    """calculate SIF using the iFLD method defined in 
    Cendrero-Mateo et al. (2019) https://doi.org/10.3390/rs11080962
    Alonso et al. (2008) https://doi.org/10.1109/LGRS.2008.2001180

    The "full_search" version finds the wavelengths one which to 
    perform the calculations itself. If processing large amounts 
    of data this can slow down the calculation significantly.
    """

    Eout_R_localmax_spect=get_local_maxima(wref_spect,sif_opts.outR_wband[0],sif_opts.outR_wband[1])
    Eout_L_localmax_spect=get_local_maxima(wref_spect,sif_opts.outL_wband[0],sif_opts.outL_wband[1])

    Eout_LRcat_localmax_spect=concat_spectra(Eout_L_localmax_spect,Eout_R_localmax_spect)

    Lout_R_localmax_spect=get_local_maxima(target_spect,sif_opts.outR_wband[0],sif_opts.outR_wband[1])
    Lout_L_localmax_spect=get_local_maxima(target_spect,sif_opts.outL_wband[0],sif_opts.outL_wband[1])
    #Lout_LRcat_localmax_spect=concat_spectra(Lout_L_localmax_spect,Lout_R_localmax_spect)

    Lout_LRcat_localmax_spect=deepcopy(Eout_LRcat_localmax_spect)
    for (n,wvl) in enumerate(Lout_LRcat_localmax_spect.wavl):
        #out=
        #print(out)
        Lout_LRcat_localmax_spect.data[n]=target_spect.closest_to_wavl(wvl)[1]

    Eout=Eout_L_localmax_spect.data[-1]    
    Eout_wvl=Eout_L_localmax_spect.wavl[-1]    

    Ein=wref_spect.min_over_band(sif_opts.in_wband)
    Ein_wvl=wref_spect.wavl_of_min_over_band(sif_opts.in_wband)

    Lout=Lout_L_localmax_spect.data[-1]
    Lin=target_spect.min_over_band(sif_opts.in_wband)

    #build interpolator for apparent reflectance 
    Rapp_spect=deepcopy(Lout_LRcat_localmax_spect)
    Rapp_spect.data/=Eout_LRcat_localmax_spect.data
    Rapp_cs=interpolate.CubicSpline(Rapp_spect.wavl,Rapp_spect.data)
    Rapp_us=interpolate.UnivariateSpline(Rapp_spect.wavl,Rapp_spect.data)

    #fit polynomial to Rapp
    coefsRapp=np.polyfit(Rapp_spect.wavl, Rapp_spect.data, 2)

        
    #fit polynomial to Eout
    coefsF=np.polyfit(Eout_LRcat_localmax_spect.wavl, Eout_LRcat_localmax_spect.data, 2)


    #plot to check E interpolation
    #(this is looking good!)
    if False:
        plt.plot(wref_spect.wavl,wref_spect.data)
        plt.plot(Eout_L_localmax_spect.wavl,Eout_L_localmax_spect.data,"o")
        plt.plot(Eout_R_localmax_spect.wavl,Eout_R_localmax_spect.data,"o")
        plt.plot(Eout_LRcat_localmax_spect.wavl,Eout_LRcat_localmax_spect.data,".")
        plt.plot(wref_spect.wavl,np.polyval(coefsF,wref_spect.wavl))
        plt.plot(Ein_wvl,np.polyval(coefs,Ein_wvl),"*")
        plt.xlim([sif_opts.outL_wband[0],sif_opts.outR_wband[1]])
        plt.show()

   
    #plot to check Rapp interpolation
    #(this shows the cubic spline looking bonkers!)
    #NOTE: DON'T USE CUBIC SPLINE - USE UNIVARIATE SPLINE
    if False:
        plt.xlim([sif_opts.outL_wband[0],sif_opts.outR_wband[1]])   
        #plt.ylim([0.0,0.2])
        plt.plot(wref_spect.wavl,target_spect.data/np.roll(wref_spect.data,0))
        plt.plot(wref_spect.wavl,Rapp_cs(wref_spect.wavl))
        plt.plot(wref_spect.wavl,Rapp_us(wref_spect.wavl),label="uni spline")
        plt.plot(wref_spect.wavl,np.polyval(coefsRapp,wref_spect.wavl))
        #plt.plot(Rapp_spect.wavl,Rapp_spect.data,"o")
        #plt.plot(wref_spect.wavl,wref_spect.data)
        #plt.plot(target_spect.wavl,target_spect.data)
        plt.legend()
        plt.show()

    
    #correction terms
    #alpha_R=Rapp_cs(Eout_wvl)/Rapp_cs(Ein_wvl)
    alpha_R=Rapp_us(Eout_wvl)/Rapp_us(Ein_wvl)
    #alpha_R=np.polyval(coefsRapp,Eout_wvl)/np.polyval(coefsRapp,Ein_wvl)
    alpha_F=alpha_R*Eout/np.polyval(coefsF,Ein_wvl)
    #alpha_R=1.
    #alpha_F=1.
    #print("* %f %f"%(alpha_R,alpha_F))
    #print("* %f %f"%(Eout_wvl,Ein_wvl))
    #print("* %f %f"%(Rapp_cs(Eout_wvl),Rapp_cs(Ein_wvl)))

    return (alpha_R*Eout*Lin-Ein*Lout)/(alpha_R*Eout-alpha_F*Ein)


def concat_spectra(spect1,spect2):
    """concatenate two spectra objects and returns
    a new spectra object.
    
    note: for speed purposes, no checking is done to 
    make sure the ordering of the wavelengths is correct.
    This is left to the user.
    """
    out=deepcopy(spect1)
    out.data=np.concatenate([spect1.data,spect2.data])
    out.wavl=np.concatenate([spect1.wavl,spect2.wavl])
    return out


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


if __name__=="__main__":

    import matplotlib.pyplot as plt

    wref=Spectra(fname="./data_in/white_reference.csv",ftype="CSV",hdrLines=0)
    test=Spectra(fname="./data_in/test_spectra.csv",ftype="CSV",hdrLines=0)

    sif_opts=sif_opts_FLD_Cendrero19_O2b
   
    if True:  
        print(sif_sFLD_full_search(test,wref,sif_opts=sif_opts)*0.01)
        print(sif_3FLD_full_search(test,wref,sif_opts=sif_opts)*0.01)
        print(sif_linReg(test,wref,sif_opts=sif_opts)*0.01)

        pseudo=sif_linReg_compute_pseudo_inverse(wref, sif_opts)
        print(sif_linReg_use_precomp_inverse(test, pseudo,sif_opts)*0.01)

        pseudo=sif_ridgeReg_compute_pseudo_inverse(wref, sif_opts=sif_opts, order=3, gamma=10000000.)
        print(sif_ridgeReg_use_precomp_inverse(test, pseudo,sif_opts=sif_opts,out_wvl=761.)*0.01)

    print(sif_iFLD_full_search(test,wref,sif_opts=sif_opts)*0.01)

    #plt.plot(wref.wavl,wref.data)
    #plt.plot(test.wavl,test.data)
    #plt.show()
    
    
    
