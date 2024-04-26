from proc_ibis import ibisBilReader

ibis_filename="./data/i200073b_mapped_utm30n.bil"
ibis_data=ibisBilReader(ibis_filename)

#white reference:
#(python and ENVI ordering are different in y)
x_wr=145
y_wr=ibis_data.lines-1137
reference=ibis_data.read_pixel(x_wr,y_wr)

with open('white_reference.csv', 'w') as f:
    for (i,data) in enumerate(reference):
        print( f"{ibis_data.wavelength[i]},{data}", file=f)

