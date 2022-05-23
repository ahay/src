import numpy as np
import os
import pylops
from curvelops import FDCT3D
import m8r
# import pyvista as pv


def byte2rsf(rsffile, inpfile, shape):
    """
    convert numpy array to rsf file
    """
    rsf = m8r.Output(rsffile)
    inp = np.array(inpfile, dtype="float32")
    rsf.put("n1", shape[2])
    rsf.put("n2", shape[1])
    rsf.put("n3", shape[0])
    rsf.put("o1", 0.0)
    rsf.put("o2", 0.0)
    rsf.put("o3", 0.0)

    rsf.write(inp)
    rsf.close()


# show data
os.system(
    """
    sfbyte < noise.rsf gainpanel=all clip=0.1 | \
    sfgrey3 flat=n frame1=80 frame2=50 frame3=50 \
    point1=0.6 point2=0.7 title= label2="X" label3="Y" unit2="km" unit3="km" | \
    sfpen
    """)

# input data
# inputfile = '../testdata/sigmoid.npz'
inp = m8r.Input("noise.rsf")
n1 = inp.int("n1")
n2 = inp.int("n2")
n3 = inp.int("n3")
d = inp.read(shape=(n3, n2, n1))

print(d.shape)
# d = np.load(inputfile)
# d = d['sigmoid']
ny, nx, nt = d.shape
dy, dx, dt = 8, 8, 0.004
y, x, t = np.arange(ny) * dy, np.arange(nx) * dx, np.arange(nt) * dt

Cop = FDCT3D(d.shape)
print(Cop.shape)


def reconstruct(data, op, perc=0.1):
    """
    Convenience function to calculate reconstruction using top
    `perc` percent of coefficients of a given operator `op`.
    """
    y = op * data.ravel()
    denoise = np.zeros_like(y)

    # Order coefficients by strength
    strong_idx = np.argsort(-np.abs(y))
    strong = np.abs(y)[strong_idx]

    # Select only top `perc`% coefficients
    strong_idx = strong_idx[:int(np.rint(len(strong_idx) * perc))]
    denoise[strong_idx] = y[strong_idx]

    data_denoise = op.inverse(denoise).reshape(data.shape)
    return np.real(data_denoise), strong


#perc = 0.08
perc = 0.1
d_dct, dct_strong = reconstruct(d, Cop, perc=perc)
dct_err = d - d_dct
print(d_dct.shape)


byte2rsf("ct_denoise.rsf", d_dct, d_dct.shape)
os.system(
    """
    sfput < ct_denoise.rsf d1=0.008 d2=0.025 d3=0.025 01=0 o2=-1.25 o3=-1.25 | \
    sfbyte gainpanel=all clip=0.1 | \
    sfgrey3 flat=n frame1=80 frame2=50 frame3=50 \
    point1=0.6 point2=0.7 title= label2="X" label3="Y" unit2="km" unit3="km" | \
    sfpen
    """)

byte2rsf("ct_denoise_err.rsf", dct_err, dct_err.shape)
os.system(
    """
    sfput < ct_denoise_err.rsf d1=0.008 d2=0.025 \
    d3=0.025 01=0 o2=-1.25 o3=-1.25 | \
    sfbyte gainpanel=all clip=0.1 | \
    sfgrey3 flat=n frame1=80 frame2=50 frame3=50 \
    point1=0.6 point2=0.7 title= label2="X" label3="Y" unit2="km" unit3="km" | \
    sfpen
    """)


# pyvista
# grid = pv.UniformGrid()
# grid.dimensions = np.array(d_dct.shape) + 1
# # Edit the spatial reference
# grid.origin = (0.0, 0.0, 0.0)  # The bottom left corner of the data set
# grid.spacing = (1, 1, 1)  # These are the cell sizes along each axis

# # Add the data values to the cell data
# grid.cell_arrays["values"] = d_dct.flatten(order="F")  # Flatten the array!
# p = pv.Plotter()
# # boxes
# bounds = [0, 50, 50, 101, 0, 126]
# clipped = grid.clip_box(bounds)
# # noninteractive
# # Set a custom position and size
# sargs = dict(
#     height=0.25,
#     vertical=True,
#     position_x=0.05,
#     position_y=0.05,
#     title_font_size=20,
#     label_font_size=16,
#     shadow=True,
#     n_labels=10,
#     italic=True,
#     fmt="%.1f",
#     font_family="arial",
#     color="black"
# )

# p = pv.Plotter()
# # p.add_mesh(grid, scalar_bar_args=sargs, cmap="seismic")
# p.add_mesh(
#     # grid,
#     clipped,
#     # color="000066",
#     # style="wireframe",
#     style="surface",
#     # show_scalar_bar=False,
#     scalar_bar_args=sargs,
#     stitle="Amplitude",
#     cmap="gray",
#     # cmap="seismic",
#     # cmap="coolwarm",
#     lighting=False,
#     # show_edges=True,
# )
# p.set_background(color="white")
# p.show_bounds(
#     xlabel="Ytrace",
#     ylabel="Xtrace",
#     zlabel="Time",
#     grid="front",
#     location="outer",
#     all_edges=True,
#     color="black"
# )
# cpos = [
#     (-296.8959883360211, 375.2608528067814, -159.3368341522591),
#     (50.0, 75.0, 100.0),
#     (0.3432324736098539, -0.3546432150194057, -0.8697238982000901),
# ]
# # p.view_vector((100, 150, 0), viewup=True)
# p.camera_position = cpos

# cpos = p.show("cmpnoise")
# # cpos = p.show("Qdome", screenshot="qdome.png")

# print(cpos)
