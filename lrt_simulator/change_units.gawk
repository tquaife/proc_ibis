#!/bin/gawk -f

BEGIN{
    #generate fluorescence file for UVSPEC 
    #converting a flat fluorescence (F) with 
    #units of mW/m2/nm to Q/s/cm2/nm
    
    #Planck constant
    h=6.626E-34
    #speed of light 
    c=3E8
}
{
    #convert wavlength to m
    w=($1*10**-9)
    #energy per photon (J/Q)
    e=h*c/w            
    #from c2 to m2
    e*=10000    
    print $1,$2*e*1000
 

}

