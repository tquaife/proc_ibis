import numpy as np
import scipy.linalg as linalg
from copy import copy
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
    """Wrapper function to preserve previous functionality
    """
    return sif_ridgeReg_compute_pseudo_inverse_RFopts(wref_spect, orderR=order, gammaR=gamma, orderF=order, gammaF=gamma,sif_opts=sif_opts)

def sif_ridgeReg_compute_pseudo_inverse_RFopts(wref_spect, orderR=1, gammaR=1., orderF=1, gammaF=1.,sif_opts=sif_opts_FLD_Cendrero19_O2a):
    """
    
    (KTK+gBTB)KT
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
    KTK=np.matmul(K.T,K)
    
    #make the block matrix for the
    #ridge regression constraints
    DR=diff_matrix(size,order=orderR)
    DF=diff_matrix(size,order=orderF)
    B=linalg.block_diag(gammaF*DF,gammaR*DR)
    BTB=np.matmul(B.T,B)
    
    return(np.matmul(np.linalg.inv(KTK+BTB),K.T))


def sif_ridgeReg_use_precomp_inverse( target_spect, pseudo_inv,sif_opts=sif_opts_FLD_Cendrero19_O2a, out_wvl=None ):

    wl_min=sif_opts.in_wband[0]
    wl_max=sif_opts.in_wband[1]

    targ=copy(target_spect)
    targ.trim(wl_min,wl_max)
    size=len(targ.data)

    if out_wvl is None:
        out_wvl=targ.wavl[int(size/2.)]

    out=np.matmul(pseudo_inv,targ.data)
    #copy SIF into targ
    targ.data=out[size:]
    
    return(targ.closest_to_wavl(out_wvl)[1])

    
    
def sif_linReg_compute_pseudo_inverse( wref_spect, sif_opts=sif_opts_FLD_Cendrero19_O2a ):

    wl_min=sif_opts.in_wband[0]
    wl_max=sif_opts.in_wband[1]

    wref=copy(wref_spect)    
    wref.trim(wl_min,wl_max)

    X=np.ones((len(wref.data),2))
    X[:,0]=np.array(wref.data)
    
    return(np.linalg.pinv(X))

def sif_linReg_use_precomp_inverse( target_spect, pseudo_inv, sif_opts=sif_opts_FLD_Cendrero19_O2a ):

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
    #wl_min=sif_opts.outL_wband[0]
    #wl_max=sif_opts.outR_wband[1]

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


if __name__=="__main__":

    import matplotlib.pyplot as plt

    wref=Spectra(fname="./data_in/white_reference.csv",ftype="CSV",hdrLines=0)
    test=Spectra(fname="./data_in/test_spectra.csv",ftype="CSV",hdrLines=0)

    sif_opts=sif_opts_FLD_Cendrero19_O2a
   
    if True:  
        print(sif_sFLD(test,wref,sif_opts=sif_opts)*0.01)
        print(sif_3FLD(test,wref,sif_opts=sif_opts)*0.01)
        print(sif_linReg(test,wref,sif_opts=sif_opts)*0.01)

        pseudo=sif_linReg_compute_pseudo_inverse(wref, sif_opts)
        print(sif_linReg_use_precomp_inverse(test, pseudo,sif_opts)*0.01)

    pseudo=sif_ridgeReg_compute_pseudo_inverse(wref, sif_opts=sif_opts, order=3, gamma=10000000.)
    print(sif_ridgeReg_use_precomp_inverse(test, pseudo,sif_opts=sif_opts,out_wvl=761.)*0.01)

    plt.plot(wref.wavl,wref.data)
    plt.plot(test.wavl,test.data)
    plt.show()
    
    
    
