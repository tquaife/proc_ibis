import numpy as np
from matplotlib import pyplot as plt
import pickle


with open('ndvi.pickle', 'rb') as f:
    image = pickle.load(f)
ndvi_data = np.ma.masked_where(image < -1000, image)
ndvi_data = ndvi_data[2000:,:]

with open('sif_o2a.pickle', 'rb') as f:
    image = pickle.load(f)
sif_o2a_data = np.ma.masked_where(image < -10, image)
sif_o2a_data = sif_o2a_data[2000:,:]

with open('sif_o2b.pickle', 'rb') as f:
    image = pickle.load(f)
sif_o2b_data = np.ma.masked_where(image < -10, image)
sif_o2b_data = sif_o2b_data[2000:,:]

#normalise:
#ndvi_data=ndvi_data/np.max(ndvi_data)
sif_o2a_data=sif_o2a_data/1.5
sif_o2b_data=sif_o2b_data/2.5

#flip axes etc for imshow
img2display=np.array([sif_o2a_data,ndvi_data,sif_o2b_data])
img2display=np.swapaxes(img2display,2,0)
img2display=np.flipud(img2display)


plt.imshow(img2display)
plt.show()



