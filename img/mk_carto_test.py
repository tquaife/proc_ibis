import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes

import cartopy
import cartopy.crs as ccrs

if __name__=="__main__":

    with open('../data_out/sif_o2a_rFLD.pickle', 'rb') as f:
        image = pickle.load(f)
    sif_o2a_data = np.ma.masked_where(image < -10, image)
    #sif_o2a_data = sif_o2a_data[2000:,:]

    (lines,samples)=np.shape(sif_o2a_data)


    #fig = plt.figure()
    #ax = fig.add_subplot(1,1,1,projection=ccrs.UTM(30))
    #ax.set_extent(extent, crs=ccrs.UTM(30))

    fig, axs = plt.subplots(nrows=1, ncols=1, sharex='col', sharey='row',
                        subplot_kw={'projection': ccrs.UTM(30)}) 


    #img_extent = (-120.67660000000001, -106.32104523100001, 13.2301484511245, 30.766899999999502)
    #img_extent = (649105.5, 649105.5+(lines*1.5), 5673408, 5673408+(samples*1.5))
    #649105.5,5673408
    #plt.imshow(sif_o2a_data,cmap="pink",vmin=-0.5,vmax=2.0,extent=img_extent,transform=ccrs.UTM(30))
    #plt.show()
