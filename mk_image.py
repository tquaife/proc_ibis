import numpy as np
from matplotlib import pyplot as plt
import pickle

with open('sif.pickle', 'rb') as f:
    image = pickle.load(f)


masked_data = np.ma.masked_where(image < -10, image)
    
plt.imshow(np.rot90(masked_data),vmin=-0.5,vmax=2.5)
plt.axis("off")
plt.colorbar(orientation="horizontal",label="SIF mW/m2/sr/um",pad=-0.0)
plt.show()
#plt.savefig("ibis_sif_example.png",bbox_inches='tight', dpi=300)
