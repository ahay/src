#!/usr/bin/env python

import mxarray as mx
import numpy as np
import xarray as xr
import sys

def test_roundtrip(path):

    # create xarray
    arr = xr.DataArray(
        np.arange(3*4*5*6, dtype=np.float32).reshape(3,4,5,6), 
        dims=("t","x","y","z"),
        coords={
            "t": [0,1,2], 
            "x": [0.0, 1.0, 2.0, 3.0], 
            "y": [0.0, 1.0, 2.0, 3.0, 4.0],
            "z": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
        }
    )
    # write to RSF
    out = path
    mx.xarray_to_rsf(arr, str(out))

    # read back
    da = mx.rsf_to_xarray(str(out))
    diff = da - arr

    assert da.shape == (3, 4, 5, 6)
    assert da.values.dtype == arr.values.dtype == np.float32
    assert da.dims == arr.dims
    assert da.coords["t"].values.tolist() == arr.coords["t"].values.tolist()
    assert da.coords["x"].values.tolist() == arr.coords["x"].values.tolist()
    assert da.coords["y"].values.tolist() == arr.coords["y"].values.tolist()
    assert da.coords["z"].values.tolist() == arr.coords["z"].values.tolist()
    assert da.values.tolist() == arr.values.tolist()
    assert diff.values.tolist() == [[[[0.0 for _ in range(6)] for _ in range(5)] for _ in range(4)] for _ in range(3)]
    
    # read as xarray Dataset
    arr = mx.rsf_to_xarrayds(path)
    
    # write to netcdf (can be loaded in paraview)
    arr.to_netcdf("test.nc")
    arr2 = xr.open_dataarray("test.nc")

    assert arr.sizes == arr2.sizes
    assert arr.coords["t"].values.tolist() == arr2.coords["t"].values.tolist()
    assert arr.coords["x"].values.tolist() == arr2.coords["x"].values.tolist()
    assert arr.coords["y"].values.tolist() == arr2.coords["y"].values.tolist()
    assert arr.coords["z"].values.tolist() == arr2.coords["z"].values.tolist()
    assert arr.data.values.tolist() == arr2.values.tolist()

    import m8r, os
    tmp = m8r.Input(file)
    binFile = tmp.string("in")
    os.remove(binFile)
    os.remove(file)
    os.remove("test.nc")

def test_cunks(path):
    arr = xr.DataArray(
        np.arange(100*100*100*5, dtype=np.float32).reshape(100,100,100,5),
        dims=("t","x","y","z"),
        coords={
            "t" : np.arange(100),
            "x" : np.arange(100),
            "y" : np.arange(100),
            "z" : np.arange(5)
        }
    )
    
    mx.xarray_to_rsf(arr, path)
    del arr
    
    import time
    start_time = time.time()
    da = mx.rsf_to_xarray(path, chunks=(2, 3, 5, 1))
    print(f"mean: {da.mean().compute().values}")
    end_time = time.time()
    print(f"read with chuks {(2, 3, 5, 1)} took {end_time - start_time} seconds")
    del da
    start_time = time.time()
    da = mx.rsf_to_xarray(path, chunks=(1, 1,1, 1))
    print(f"mean: {da.mean().compute().values}")
    end_time = time.time()
    print(f"read with chuks {(1, 1, 1, 1)} took {end_time - start_time} seconds")
    del da
    start_time = time.time()
    da = mx.rsf_to_xarray(path, chunks="auto")
    print(f"mean: {da.mean().compute().values}")
    end_time = time.time()
    print(f"read with chuks auto took {end_time - start_time} seconds")
    del da
    
    import os
    os.remove(path)
    

if __name__ == "__main__":
    file = "tmp.rsf"
    test_roundtrip(path=file)
    # test_cunks(path="./tmp_large.rsf")
    sys.stdout.write("All tests passed.\n")
## test chunks
####################################################################################
############# output for data with shape (100, 100, 100, 5) ########################
####################################################################################
# mean: 2499999.5
# read with chuks (2, 3, 5, 1) took 21.33984923362732 seconds
# mean: 2499999.5
# read with chuks (1, 1, 1, 1) took 559.3831532001495 seconds
# mean: 2499999.5
# read with chuks auto took 0.01817178726196289 seconds
####################################################################################