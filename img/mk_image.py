import pickle
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

with open('../data_out/sif_o2a_ridgeReg2.pickle', 'rb') as f:
    image = pickle.load(f)
sif_o2a_data = np.ma.masked_where(image < -10, image)
sif_o2a_data = sif_o2a_data[2000:,:]

with open('../data_out/sif_o2b_ridgeReg2.pickle', 'rb') as f:
    image = pickle.load(f)
sif_o2b_data = np.ma.masked_where(image < -10, image)
sif_o2b_data = sif_o2b_data[2000:,:]

cmap="pink"
#cmap="hot"

fig, axes = plt.subplots(nrows=2, ncols=1)

#O2B ========================
vmin=-0.2
vmax=1.0
im0=axes[0].imshow(np.rot90(sif_o2b_data), vmin=vmin, vmax=vmax, cmap=cmap)
box = axes[0].get_position()
#print(box)
#box.y0 = box.y0 - 0.2
#box.y1 = box.y1 - 0.2
#axes[0].set_position(box)
axes[0].axis("off")

divider = make_axes_locatable(axes[0])
cax = divider.append_axes('bottom', size='10%', pad=0.05)
fig.colorbar(im0, cax=cax, orientation="horizontal", label="O2b SIF mW/m2/sr/um",pad=0.01, cmap=cmap)


#O2A ========================
vmin=-0.4
vmax=2.0
im1=axes[1].imshow(np.rot90(sif_o2a_data), vmin=vmin, vmax=vmax, cmap=cmap)
axes[1].axis("off")
fig.subplots_adjust(wspace=0, hspace=-0.4)
divider = make_axes_locatable(axes[1])
cax = divider.append_axes('bottom', size='10%', pad=0.05)
fig.colorbar(im1, cax=cax, orientation="horizontal", label="O2a SIF mW/m2/sr/um",pad=0.01, cmap=cmap)


#axes[0].text(0.5, 0.8, 'SIF O2b', horizontalalignment='center', verticalalignment='center', transform=axes[0].transAxes)
#axes[1].text(0.5, 0.8, 'SIF O2a', horizontalalignment='center', verticalalignment='center', transform=axes[1].transAxes)

#cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.065])
#fig.colorbar(im, orientation="horizontal", label="SIF mW/m2/sr/um",pad=0.01, cmap=cmap, cax=cbar_ax)

#plt.show()
fig.tight_layout()

plt.savefig("ibis_sif_ridgeReg2_o2a_o2b.png",bbox_inches='tight', dpi=300)


