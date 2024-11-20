import sys
sys.path.append("../")

import numpy as np
from matplotlib import pyplot as plt
from spectra_tools.spectra import Spectra

def convert_units_Qcm2_to_Wm2(wvl,datum):
    """
    input should be in Q/s/cm2/nm
    n.b. Q is *not* moles, it is the number of photons
    
    output is in W/m2/nm
    """

    #Planck constant
    h=6.626E-34
    #speed of light 
    c=2.99792E8
    #convert wavlength to m
    w=(wvl*10**-9)
    #energy per photon (J/Q)
    e=h*c/w            
    #from c2 to m2
    e*=10000
    #convert data to W/m2
    return datum*e


class libRadSpectra(Spectra):

    def __init__(self, fname=None,ftype="TXT",wavlCol=0,dataCol=1,hdrLines=0):
        super().__init__(fname=fname,ftype=ftype,wavlCol=wavlCol,dataCol=dataCol,hdrLines=hdrLines)

    def convert_units_Qcm2_to_Wm2(self):
        """
        Applies convert_units_Qcm2_to_Wm2 to the entire spectrum

        input should be in Q/s/cm2/nm
        n.b. Q is *not* moles, it is the number of photons
        
        output is in W/m2/nm
       
        do not apply this method multiple times to the data
        could really do with some checking of units in the 
        base class.
        """
        
        for (i,w) in enumerate(self.wavl):
            self.data[i]=convert_units_Qcm2_to_Wm2(w,self.data[i])
    
         
    def resample_to_ibis(self,ibis_wavls,band_width=0.35): #original: band_width=0.11
        """resample a libradtran spectra to IBIS bands 
        
        ibis_wavls is an array of band centers
        band_width is the wifht of the ibis bands
        
        n.b. wavelength units should match with
        whatever has been produced by libradtran
        (typically nm).
        
        the difference between band centres is 0.11nm but
        averaging over a band width of 0.11nm gives spectra
        that are not nearly as smooth as the IBIS.
        Using 0.35nm gives a more realistic result.
        """
        data_new=[]
        for w in ibis_wavls:
             tmp=self.data[np.abs(self.wavl-w)<=band_width/2.]
             data_new.append(tmp.mean())

        self.wavl=np.array(ibis_wavls)
        self.data=np.array(data_new)




if __name__=="__main__":

    fn="lrt_io/uvspec_fluorescence_z1000m.out"
    #s=libRadSpectra(fn,ftype="TXT",dataCol=2,hdrLines=0)
    s=libRadSpectra(fn)
    print(s.ftype)
    s.convert_units_Qcm2_to_Wm2()
    #plt.plot(s.wavl,s.data)    
    ibis_wavls=np.genfromtxt("ibis_wavelengths.txt")    
    s.resample_to_ibis(ibis_wavls)

    fn="lrt_io/uvspec_fluorescence_z1000m_no_fluor.out"
    n=libRadSpectra(fn,ftype="TXT",dataCol=2,hdrLines=0)
    n.convert_units_Qcm2_to_Wm2()
    #plt.plot(s.wavl,s.data)    
    ibis_wavls=np.genfromtxt("ibis_wavelengths.txt")    
    n.resample_to_ibis(ibis_wavls)


    plt.plot(s.wavl,s.data,label="w/ SIF")
    plt.plot(n.wavl,n.data,label="no SIF")
    #plt.plot(n.wavl,s.data-n.data,label="difference")
    plt.legend()
    plt.show()
    
    
