##   xarray interface for RSF files
##
##   Copyright (C) 2025 University of Texas at Austin
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import m8r
import numpy as np

try:
    import xarray as xr
except ImportError:
    xr = None
    
try:
    import dask.array as da
except ImportError:
    da = None

def rsf_to_xarray(path, chunks = "auto"):
    """Convert an RSF file to an xarray DataArray with optimal Dask lazy loading."""
    
    if xr is None:
        raise ImportError("xarray is required.")
    
    f = m8r.Input(path)
    
    shape = f.shape()
    ndim = len(shape)
    
    dims = []
    # units = []
    coords = {}
    for axis in range(1, ndim+1):
        n = f.int(f"n{axis}")
        o = f.float(f"o{axis}")
        d = f.float(f"d{axis}")
        label = f.string(f"label{axis}")
        if label == None:
            label = f"dim{axis}"
        # unit = f.string(f"unit{axis}")
        
        coords[label] = np.arange(n) * d + o
        dims.append(label)
        # units.append(unit)
        
    # data = np.asarray(f)
    binFile = f.string("in")
    dtype = f.string("type")
    if dtype is None or dtype == "float":
        dtype = np.float32
    elif dtype == "int":
        dtype = np.int32
    else:
        raise ValueError(f"Unsupported data type: {dtype}")
    
    mm = np.memmap(
        binFile,
        dtype=dtype,
        mode='r',
        shape=shape
    )

    if da is not None:
        data = da.from_array(mm, chunks=chunks)
    else:
        data = np.asarray(mm)
        

    # covert c order to python order
    # data = data.reshape(shape_py)
    data = data.transpose(*reversed(range(data.ndim)))

    ds = xr.DataArray(
        data, 
        dims=dims, 
        coords=coords
    )
    
    return ds

def xarray_to_rsf(ds, outpath):
    """Convert an xarray Dataset to an RSF file."""
    
    if xr is None:
        raise ImportError("xarray is required.")
    
    if not isinstance(ds, xr.DataArray):
        raise ValueError("Input must be an xarray DataArray.")
    
    data = ds.values
    data = data.transpose(*reversed(range(data.ndim)))
    
    dims = ds.dims

    
    out = m8r.Output(outpath)
    
    for i, dim in enumerate(dims, start=1):
            
        coord = ds.coords[dim].values
        if len(coord) > 1:
            d = coord[1] - coord[0]
        else:
            d = np.float32(1.)
            
        o = coord[0]
        n = len(coord)
        
        out.put(f"n{i}", n)
        out.put(f"o{i}", o)
        out.put(f"d{i}", d)
        out.put(f"label{i}", str(dim))
        
    out.write(data.astype(np.float32))
    out.close()
    
def rsf_to_xarrayds(path, chunks = "auto"):
    """Convert an RSF file to an xarray Dataset."""
    
    da = rsf_to_xarray(path, chunks=chunks)
    ds = da.to_dataset(name="data")
    
    return ds