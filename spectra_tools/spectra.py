#!/usr/bin/env python

import re
from copy import copy
import numpy as np


class UnknownFileType(Exception):
  """Exception class for unknown filetypes
  """
  pass

class Spectra(object):

  def __init__(self,fname=None,ftype="SVC",wavlCol=0,dataCol=1,hdrLines=1):
    
    self.data=np.array([])
    self.wavl=np.array([])
    self.ftype=ftype
    
    if fname!=None:
      self.loadSpectra(fname,wavlCol,dataCol,hdrLines)
      
  def loadSpectra(self,fname,wavlCol=0,dataCol=1,hdrLines=1):    
    """Load in the spectra from a given file using a
    method appropriate to the type of file.
    
    Arguments:
    
    fname    - valid filename containing spectra
    
    Supported formats:
    
    SVC      - SCV .sig ascii file
    CSV      - standard ascii comma seperated values
    """
    
    with open(fname) as f:
      if self.ftype=="SVC":
        self.loadSVCSig(f)   
      elif self.ftype=="CSV":
        self.loadCSV(f,wavlCol,dataCol,hdrLines) 
      else:
        raise UnknownFileType(self.ftype)
      f.close()
    
  def loadCSV(self,f,wavlCol=0,dataCol=1,hdrLines=1):
    """Read in data from a standard CSV file
    
    Arguments:
    
    f        - File object
    wavlCol  - column containing wavelengths
    dataCol  - column containing dataectance data
    hdrLines - Number of lines to skip at start of file
    """
    
    tmp=np.loadtxt(f,delimiter=",",skiprows=hdrLines,usecols=(wavlCol,dataCol))
    self.wavl=tmp[:,0]
    self.data=tmp[:,1]
    
  def loadSVCSig(self,f):  
    """Read in data from an SVC .sig ascii file

    Arguments:
    
    f        - File object
    """
    getData=False

    for line in f:
      if getData:
        dataTmp=np.append(self.data, float(line.split()[3])/100.)
        wavlTmp=np.append(self.wavl, float(line.split()[0]))
        self.data=copy(dataTmp)
        self.wavl=copy(wavlTmp)
      if re.match('data=',line):
        getData=True

  def avg_over_band(self, wband):
    """calculates the mean data over the waveband 
    specified by wmin and wmax 
    """
    wmin=wband[0]
    wmax=wband[1]
    return np.mean(self.data[(self.wavl>=wmin)&(self.wavl<=wmax)])
     
  def max_over_band(self, wband):
    """calculates the maximum data value over the waveband 
    specified by wmin and wmax 
    """
    wmin=wband[0]
    wmax=wband[1]
    return np.max(self.data[(self.wavl>=wmin)&(self.wavl<=wmax)])
     
  def min_over_band(self, wband):
    """calculates the minimum data value over the waveband 
    specified by wmin and wmax 
    """
    wmin=wband[0]
    wmax=wband[1]
    return np.min(self.data[(self.wavl>=wmin)&(self.wavl<=wmax)])

            
  def interpolate(self,resltn=0.1):
    """Interpolate spectra to the given resolution.
    Overwites existing data.
    
    Arguments:
    
    resltn    - resolution of the interpolation
    """      
    
    #find the starting and ending wavelengths    
    begWavl=np.ceil(self.wavl[0]/resltn)*resltn
    endWavl=np.floor(self.wavl[-1]/resltn)*resltn

    #print self.wavl[0], begWavl
    #print self.wavl[-1], endWavl

    #generate new wavelength and relfectance arrays
    wavlTmp=np.arange(begWavl, endWavl+resltn, resltn)
    dataTmp=np.zeros(np.shape(wavlTmp))

    #perfrom a linear interpolation:
    m=0
    for (n,wavl) in enumerate(wavlTmp):

       while self.wavl[m]<wavl and wavl != wavlTmp[-1]:
         m+=1
         
       if self.wavl[m]==wavl:
         dataTmp[n]=self.data[m]
       else:
         w1=self.wavl[m-1]
         w2=self.wavl[m]
         r1=self.data[m-1]
         r2=self.data[m]
         f=(w2-wavl)/(w2-w1)
         dataTmp[n]=r1*f+r2*(1-f)
    
    #copy in interploated data             
    self.wavl=copy(wavlTmp)
    self.data=copy(dataTmp)  
   
  def trim(self,wlmin,wlmax):
    """Trim the spectra so it is between
    two specified wavelengths. Destroys
    the original data.

    Arguments:
    
    wlmin    - the lowest wavelength of the new spectra
    wlmax    - the highest wavelength of the new spectra
    """
    
    dataTmp=np.array([])
    wavlTmp=np.array([])

    for (n,wavl) in enumerate(self.wavl):
       if (wavl>=wlmin-1e-09) and (wavl<=wlmax+1e-09): 
         wavlTmp2=np.append(wavlTmp,self.wavl[n])
         dataTmp2=np.append(dataTmp,self.data[n])
         wavlTmp=copy(wavlTmp2)
         dataTmp=copy(dataTmp2)
    #copy over trimmed data             
    self.wavl=copy(wavlTmp)
    self.data=copy(dataTmp)  

    
def convolve(s1orig,s2orig,resln=1.0,s2norm=True):
  """Convolve one spectra with another, for example
  to apply a band pass, or a spectral response function.
  
  Arguments:
  
  s1      - a spectra object
  s2      - a spectra object
  resln   - the spectral resolution to use
  s2norm  - if True normalise the second spectra 
            (e.g. to apply a spectra response function).
  """  
    
  #make copies so as not to alter
  #original data
  s1=copy(s1orig)
  s2=copy(s2orig)

  #interpolate to common resolution
  s1.interpolate(resln)
  s2.interpolate(resln)

    
  #trim spectra to encompass the exclusive
  #range of the two
  wlmin=np.max([s1.wavl[0],s2.wavl[0]])  
  wlmax=np.min([s1.wavl[-1],s2.wavl[-1]])  
  s1.trim(wlmin,wlmax)
  s2.trim(wlmin,wlmax)

  #convolve and normailse if required
  norm=1.0
  if s2norm:
    norm=s2.data.sum()
  return np.dot(s1.data,s2.data)/norm
    

def sentinel2(s,mission="a"):

  srf=[]
  sen2=spectra()
  
  if mission=="a":  
    fname="S2a_SRF.csv"
  else:
    fname="S2b_SRF.csv"
  
  for n in xrange(1,14):
    srf=spectra(fname=fname,ftype="CSV",dataCol=n)
    dataTmp=np.append(sen2.data,convolve(s,srf,0.1))
    sen2.data=copy(dataTmp)
    m=np.argmax(srf.data)
    wavlTmp=np.append(sen2.wavl,srf.wavl[m])
    sen2.wavl=copy(wavlTmp)
  
  return sen2
  
    
if __name__=="__main__":

  from matplotlib import pyplot as plt

  doTest1=True
  doTest2=False

  if doTest1:
    #test simulation of S2 bands
    svc=spectra(fname="HRPDA.053017.0065_moc.sig")
    sen2a=sentinel2(svc)
    sen2b=sentinel2(svc,mission="b")

    plt.plot(svc.wavl,svc.data)
    plt.plot(sen2a.wavl,sen2a.data,'o-')
    plt.plot(sen2b.wavl,sen2b.data,'o-')
    plt.xlabel('wavelength (nm)')
    plt.ylabel('relfectance (-)')
    plt.show()
    
  if doTest2:
    #test interpolation and trim routines
    s=spectra(fname="HRPDA.053017.0065_moc.sig")
    plt.plot(s.wavl,s.data)
    s.interpolate(50)
    plt.plot(s.wavl,s.data,'o')
    s.interpolate(100)
    s.trim(1200,2000)
    plt.plot(s.wavl,s.data,'--')
    plt.show()  

