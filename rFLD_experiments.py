from proc_ibis import *


sif_opts_rFLD_custom_O2a=sif_opts_base()
sif_opts_rFLD_custom_O2b=sif_opts_base()
#sif_opts_rFLD_custom_O2a.in_wband=(759.5, 761.5) 
#sif_opts_rFLD_custom_O2b.in_wband=(685.5, 688) 
#sif_opts_rFLD_custom_O2a.in_wband=(759.5, 765) 
#sif_opts_rFLD_custom_O2b.in_wband=(686, 697) 


if __name__=="__main__":

    #open ibis data 
    ibis_filename="./data_in/i200073b_mapped_utm30n.bil"
    ibis_data=ibisBilReader(ibis_filename)

    #read in reference data
    wref=Spectra(fname="./data_in/white_reference.csv",ftype="CSV",hdrLines=0)

    #rFLD O2a 
    out_fn='./data_out/sif_o2a_rFLD.pickle'
    
    #sif_opts=sif_opts_rFLD_custom_O2a
    sif_opts=sif_opts_FLD_Cendrero19_O2a

    #ibis_sif_rFLD(wref,ibis_data,3,1000,2,100000,sif_opts,out_fn,761.0) 
    #ibis_sif_rFLD(wref,ibis_data,2,10000,2,10000,sif_opts,out_fn,761.0)     
    ibis_sif_rFLD(wref,ibis_data,1,1,1,1,sif_opts,out_fn,761.0)     

    #rFLD O2b
    out_fn='./data_out/sif_o2b_rFLD.pickle'
    sif_opts=sif_opts_rFLD_custom_O2b
    sif_opts=sif_opts_FLD_Cendrero19_O2b
    #ibis_sif_rFLD(wref,ibis_data,2,0.001,2,1000,sif_opts,out_fn,686.7) 




    #linear FLD O2a
    out_fn='./data_out/sif_o2a_linFLD.pickle'
    sif_opts=sif_opts_rFLD_custom_O2a
    #sif_opts=sif_opts_FLD_Cendrero19_O2a
    #ibis_sif_linearFLD(wref,ibis_data,sif_opts,out_fn) 

    #linear FLD O2b
    out_fn='./data_out/sif_o2b_linFLD.pickle'
    sif_opts=sif_opts_rFLD_custom_O2b
    #sif_opts=sif_opts_FLD_Cendrero19_O2b
    #ibis_sif_linearFLD(wref,ibis_data,sif_opts,out_fn) 



    #rFLD O2a testing wavebands and gamma
    out_fn='./data_out/sif_o2a_rFLD.pickle'
    sif_opts=sif_opts_rFLD_custom_O2a
    #ibis_sif_rFLD(wref,ibis_data,3,1000,2,100000,sif_opts,out_fn,761.0) 
    #ibis_sif_rFLD(wref,ibis_data,3,0.01,2,1,sif_opts,out_fn,761.0)     

    #rFLD O2b
    out_fn='./data_out/sif_o2b_rFLD.pickle'
    #sif_opts=sif_opts_rFLD_custom_O2b
    #ibis_sif_rFLD(wref,ibis_data,2,0.001,2,1000,sif_opts,out_fn,686.7) 



    #rFLD O2a - custom waveband range:
    #out_fn='./data_out/sif_o2a_rFLD_custom.pickle'
    #sif_opts=sif_opts_rFLD_custom_O2a
    #ibis_sif_rFLD(wref,ibis_data,3,10E6,2,10E6,sif_opts,out_fn,761.0) 
    

    #rFLD O2a
    #out_fn='./data_out/sif_o2a_rFLD_32.pickle'
    #sif_opts=sif_opts_FLD_Cendrero19_O2a
    #ibis_sif_rFLD(wref,ibis_data,3,10E6,2,10E6,sif_opts,out_fn,761.0) 

    #rFLD O2b
    #out_fn='./data_out/sif_o2b_rFLD_32.pickle'
    #sif_opts=sif_opts_FLD_Cendrero19_O2b
    #ibis_sif_rFLD(wref,ibis_data,3,10E6,2,10E6,sif_opts,out_fn,686.7) 

