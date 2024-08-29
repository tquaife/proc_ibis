from proc_ibis import *

if __name__=="__main__":

    #open ibis data 
    ibis_filename="./data_in/i200073b_mapped_utm30n.bil"
    ibis_data=ibisBilReader(ibis_filename)

    #read in reference data
    wref=Spectra(fname="./data_in/white_reference.csv",ftype="CSV",hdrLines=0)

    #rFLD O2a
    out_fn='./data_out/sif_o2a_rFLD_32.pickle'
    sif_opts=sif_opts_FLD_Cendrero19_O2a
    ibis_sif_rFLD(wref,ibis_data,3,10E6,2,10E6,sif_opts,out_fn,761.0) 

    #rFLD O2b
    out_fn='./data_out/sif_o2b_rFLD_32.pickle'
    sif_opts=sif_opts_FLD_Cendrero19_O2b
    ibis_sif_rFLD(wref,ibis_data,3,10E6,2,10E6,sif_opts,out_fn,686.7) 

