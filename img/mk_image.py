import pickle
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def mk_img(name):

    with open('../data_out/sif_o2a_%s.pickle'%name, 'rb') as f:
        image = pickle.load(f)
    sif_o2a_data = np.ma.masked_where(image < -10, image)
    #sif_o2a_data = sif_o2a_data[2000:,:]
    sif_o2a_data = sif_o2a_data[:500,:]

    with open('../data_out/sif_o2a_%s.pickle'%name, 'rb') as f:
        image = pickle.load(f)
    sif_o2b_data = np.ma.masked_where(image < -10, image)
    #sif_o2b_data = sif_o2b_data[2000:,:]
    sif_o2b_data = sif_o2b_data[:500,:]

    cmap="pink"
    #cmap="hot"

    fig, axes = plt.subplots(nrows=2, ncols=1)

    #O2B ========================
    vmin=-0.2
    vmax=1.5
    #vmax=3.5
    im0=axes[0].imshow(np.rot90(sif_o2b_data), vmin=vmin, vmax=vmax, cmap=cmap)
    axes[0].axis("off")
    divider = make_axes_locatable(axes[0])
    cax = divider.append_axes('bottom', size='10%', pad=0.05)
    fig.colorbar(im0, cax=cax, orientation="horizontal", label="O2b SIF mW/m2/sr/um",pad=0.01, cmap=cmap)

    #O2A ========================
    vmin=-0.4
    vmax=2.0
    #vmax=4.0
    im1=axes[1].imshow(np.rot90(sif_o2a_data), vmin=vmin, vmax=vmax, cmap=cmap)
    axes[1].axis("off")
    fig.subplots_adjust(wspace=0, hspace=-0.4)
    divider = make_axes_locatable(axes[1])
    cax = divider.append_axes('bottom', size='10%', pad=0.05)
    fig.colorbar(im1, cax=cax, orientation="horizontal", label="O2a SIF mW/m2/sr/um",pad=0.01, cmap=cmap)

    #save figure:
    fig.tight_layout()
    plt.savefig("ibis_sif_%s_o2a_o2b.png"%name,bbox_inches='tight', dpi=300)


if __name__=="__main__":

    proc_list=[]
    #proc_list.append("sFLD_fast")
    #proc_list.append("sFLD")
    proc_list.append("iFLD")
    #proc_list.append("3FLD_fast")
    #proc_list.append("3FLD")
    #proc_list.append("ridgeReg")
    #proc_list.append("linReg")
    
    for name in proc_list:
        mk_img(name)



