from proc_ibis import *
"""Process IBIS data to retrieve SIF.

This is the simplest possible example of using the rFLD
with a white reference and retrieving O2a SIF at 761nm.
Results are written to a pickle file.
"""
#open ibis data 
ibis_filename="./data_in/i200073b_mapped_utm30n.bil"
ibis_data=ibisBilReader(ibis_filename)

#read in reference data
wref=Spectra(fname="./data_in/white_reference.csv",ftype="CSV",hdrLines=0)

#rFLD O2a
out_fn='./data_out/sif_o2a_rFLD.pickle'
sif_opts=sif_opts_FLD_Cendrero19_O2a
ibis_sif_sFLD(wref,ibis_data,2,10E6,2,10E6,sif_opts,out_fn,761.0) 



