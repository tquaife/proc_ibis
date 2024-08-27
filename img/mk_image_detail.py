import pickle
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable



def mk_img(band):

    #band="o2a"
    im_xmin=2200
    im_xmax=2750
    im_ymin=20
    im_ymax=420
    cmap="pink"
    
    if band =="o2a":
        vmin=-0.4
        vmax=2.0
    elif band =="o2b":
        vmin=-0.2   
        vmax=1.5

    with open('../data_out/sif_%s_ridgeReg.pickle'%band, 'rb') as f:
        image = pickle.load(f)
    sif_rFLD_data = np.ma.masked_where(image < -10, image)
    sif_rFLD_data = np.rot90(sif_rFLD_data[im_xmin:im_xmax,im_ymin:im_ymax])

    with open('../data_out/sif_%s_sFLD.pickle'%band, 'rb') as f:
        image = pickle.load(f)
    sif_sFLD_data = np.ma.masked_where(image < -10, image)
    sif_sFLD_data = np.rot90(sif_sFLD_data[im_xmin:im_xmax,im_ymin:im_ymax])

    with open('../data_out/sif_%s_3FLD.pickle'%band, 'rb') as f:
        image = pickle.load(f)
    sif_3FLD_data = np.ma.masked_where(image < -10, image)
    sif_3FLD_data = np.rot90(sif_3FLD_data[im_xmin:im_xmax,im_ymin:im_ymax])

    with open('../data_out/sif_%s_iFLD.pickle'%band, 'rb') as f:
        image = pickle.load(f)
    sif_iFLD_data = np.ma.masked_where(image < -10, image)
    sif_iFLD_data = np.rot90(sif_iFLD_data[im_xmin:im_xmax,im_ymin:im_ymax])

    fig, axes = plt.subplots(nrows=2, ncols=2)

    axes[0,0].imshow(sif_sFLD_data, vmin=vmin, vmax=vmax, cmap=cmap)
    axes[0,0].axis("off")
    axes[0,0].set_title("sFLD")

    axes[0,1].imshow(sif_3FLD_data, vmin=vmin, vmax=vmax, cmap=cmap)
    axes[0,1].axis("off")
    axes[0,1].set_title("3FLD")

    axes[1,0].imshow(sif_iFLD_data, vmin=vmin, vmax=vmax, cmap=cmap)
    axes[1,0].axis("off")
    axes[1,0].set_title("iFLD")

    im=axes[1,1].imshow(sif_rFLD_data, vmin=vmin, vmax=vmax, cmap=cmap)
    axes[1,1].axis("off")
    axes[1,1].set_title("rFLD")

    fig.subplots_adjust(bottom=0.2)
    #fig.tight_layout()
    cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.065])
    fig.colorbar(im, orientation="horizontal", label="SIF mW/m2/sr/um",pad=0.01, cmap=cmap, cax=cbar_ax)
    
    #plt.show()
    plt.savefig("ibis_sif_detail_%s.png"%band,bbox_inches='tight', dpi=300)
   
    
if __name__=="__main__":
    mk_img("o2a")    
    mk_img("o2b")    
    
    
