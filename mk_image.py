import numpy as np
from matplotlib import pyplot as plt
import pickle

with open('sif_o2a.pickle', 'rb') as f:
    image = pickle.load(f)
sif_o2a_data = np.ma.masked_where(image < -10, image)
sif_o2a_data = sif_o2a_data[2000:,:]

with open('sif_o2b.pickle', 'rb') as f:
    image = pickle.load(f)
sif_o2b_data = np.ma.masked_where(image < -10, image)
sif_o2b_data = sif_o2b_data[2000:,:]



vmin=-0.5
vmax=2.0


fig, axes = plt.subplots(nrows=2, ncols=1)
axes[0].imshow(np.rot90(sif_o2b_data), vmin=vmin, vmax=vmax)
axes[0].axis("off")
im=axes[1].imshow(np.rot90(sif_o2a_data), vmin=vmin, vmax=vmax)
axes[1].axis("off")
fig.subplots_adjust(wspace=0, hspace=-0.4)

axes[0].text(0.5, 0.8, 'SIF O2b', horizontalalignment='center', verticalalignment='center', transform=axes[0].transAxes)
axes[1].text(0.5, 0.8, 'SIF O2a', horizontalalignment='center', verticalalignment='center', transform=axes[1].transAxes)

cbar_ax = fig.add_axes([0.15, 0.1, 0.7, 0.065])
fig.colorbar(im, orientation="horizontal", label="SIF mW/m2/sr/um",pad=0.01, cax=cbar_ax)

#plt.show()
plt.savefig("ibis_sif_sFLD_o2a_o2b.png",bbox_inches='tight', dpi=300)

#plt.subplot(2, 1, 1)
#plt.imshow(np.rot90(sif_o2b_data), vmin=vmin, vmax=vmax)
#plt.axis("off")
#plt.subplot(2, 1, 2)
#plt.imshow(np.rot90(sif_o2a_data), vmin=vmin, vmax=vmax)
#plt.axis("off")
#plt.colorbar(orientation="horizontal",label="SIF mW/m2/sr/um",pad=0.01)
plt.show()

#plt.imshow(np.rot90(sif_o2a_data),vmin=-,vmax=2.5)
#plt.axis("off")
#plt.colorbar(orientation="horizontal",label="SIF mW/m2/sr/um",pad=0.01)
#plt.show()

