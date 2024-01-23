import numpy as np
import os
import pylops
from curvelops import FDCT2D
import m8r


def byte2rsf(rsffile, inpfile, shape):
    """
    convert numpy array to rsf file
    """
    rsf = m8r.Output(rsffile)
    inp = np.array(inpfile, dtype="float32")
    rsf.put("n1", shape[1])
    rsf.put("n2", shape[0])
    rsf.put("o1", 0.0)
    rsf.put("o2", 0.0)

    rsf.write(inp)
    rsf.close()


# show data
os.system("""
    sfgrey < data.rsf title= | sfpen
    """)

inp = m8r.Input("data.rsf")
n1 = inp.int("n1")
n2 = inp.int("n2")
d = inp.read(shape=(n2, n1))
# print(d.shape)

nx, nt = d.shape
dx, dt = 0.04238, 0.004
x, t = np.arange(nx) * dx, np.arange(nt) * dt

Cop = FDCT2D(d.shape)
Wop = pylops.signalprocessing.DWT2D(d.shape, wavelet="db9", level=2)
# print(Cop.shape)


def reconstruct(data, op, perc=0.1, inv=False):
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
    if inv:
        data_denoise = op.inverse(denoise).reshape(data.shape)
    else:
        data_denoise = (op.H * denoise).reshape(data.shape)
    return np.real(data_denoise), strong


#perc = 0.08
cperc = 0.06
d_dct, dct_strong = reconstruct(d, Cop, perc=cperc, inv=True)
# dct_err = d - d_dct
# print(d_dct.shape)

byte2rsf("ct_denoise.rsf", d_dct, d_dct.shape)
os.system("""
    sfput < ct_denoise.rsf d1=0.008 d2=0.04238 o1=0 o2=0  | \
    sfgrey title= clip=3.10373 font=2 titlefat=2 labelfat=3 gridfat=3 | sfpen
    """)