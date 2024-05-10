import sys
sys.path.append("../")

import numpy as np
from matplotlib import pyplot as plt
from spectra_tools.spectra import Spectra

class libRadSpectra(Spectra):

    def convert_units_Qcm2_to_Wm2(self):
        """
        input should be in Q/s/cm2/nm
        n.b. Q is *not* moles, it is the number of photons
        
        output is in W/m2/nm
       
        do not apply this method multiple times to the data
        could really do with some checking of units in the 
        base class.
        """
                
        #Planck constant
        h=6.626E-34
        #speed of light 
        c=3E8

        for (i,w) in enumerate(self.wavl):
            #convert wavlength to m
            w=(w*10**-9)
            #energy per photon (J/Q)
            e=h*c/w            
            #from c2 to m2
            e*=10000
            #convert data to W/m2
            self.data[i]=self.data[i]*e
         
    def resample_to_ibis(self,ibis_wavls,band_width=0.11):
        """resample a libradtran spectra to IBIS bands 
        
        ibis_wavls is an array of band centers
        band_width is the wifht of the ibis bands
        
        n.b. wavelength units should match with
        whatever has been produced by libradtran
        (typically nm).
        """
        data_new=[]
        for w in ibis_wavls:
             tmp=self.data[np.abs(self.wavl-w)<=band_width/2.]
             data_new.append(tmp.mean())

        self.wavl=np.array(ibis_wavls)
        self.data=np.array(data_new)




if __name__=="__main__":

    fn="UVSPEC_FLUORESCENCE.OUT"
    s=libRadSpectra(fn,ftype="TXT",dataCol=2,hdrLines=0)
    s.convert_units_Qcm2_to_Wm2()
    
    plt.plot(s.wavl,s.data)
    
    ibis_wavls=np.genfromtxt("ibis_wavelengths.txt")    
    s.resample_to_ibis(ibis_wavls)

    plt.plot(s.wavl,s.data)

    plt.show()
    
    
