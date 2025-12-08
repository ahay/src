#!/usr/bin/env python

import unittest
import numpy as np
import xarray as xr
import m8r
import os
import sys
import time
import re

class TestMXArray(unittest.TestCase):

    def setUp(self):
        """Runs before every test."""
        self.rsf_path = "tmp_test.rsf"
        self.nc_path = "tmp_test.nc"
        self.files_to_clean = [self.rsf_path, self.nc_path]

    def tearDown(self):
        """Runs after every test. Clean up RSF headers AND binaries."""
        
        if os.path.exists(self.rsf_path):
            try:
                with open(self.rsf_path, 'r', errors='ignore') as f:
                    content = f.read()
                    # Regex to find in="/path/to/binary" or in=binary
                    match = re.search(r'in="?([^"\s]+)"?', content)
                    if match:
                        bin_file = match.group(1)
                        if bin_file not in ['stdin', 'stdout'] and os.path.exists(bin_file):
                            os.remove(bin_file)
            except Exception as e:
                print(f"Warning: Could not clean up binary for {self.rsf_path}: {e}")

        for f in self.files_to_clean:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError:
                    pass
                    
        if hasattr(self, 'generated_temps'):
            for temp in self.generated_temps:
                if os.path.exists(temp): os.remove(temp)

    def test_roundtrip(self):
        """Test writing xarray to RSF and reading it back."""
        # Create xarray
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

        # Write to RSF
        mx.xarray_to_rsf(arr, self.rsf_path)

        # Read back as DataArray
        da = mx.rsf_to_xarray(self.rsf_path)
        diff = da - arr

        # Assertions
        self.assertEqual(da.shape, (3, 4, 5, 6))
        self.assertEqual(da.values.dtype, np.float32)
        self.assertEqual(da.dims, arr.dims)
        np.testing.assert_array_equal(da.coords["t"].values, arr.coords["t"].values)
        np.testing.assert_array_equal(da.coords["x"].values, arr.coords["x"].values)
        np.testing.assert_array_equal(da.coords["y"].values, arr.coords["y"].values)
        np.testing.assert_array_equal(da.coords["z"].values, arr.coords["z"].values)
        np.testing.assert_array_equal(da.values, arr.values)
        
        # Verify difference is zero
        expected_diff = np.zeros((3,4,5,6), dtype=np.float32)
        np.testing.assert_array_equal(diff.values, expected_diff)

        # Test Dataset conversion and NetCDF export
        ds = mx.rsf_to_xarrayds(self.rsf_path)
        ds.to_netcdf(self.nc_path)
        
        arr2 = xr.open_dataarray(self.nc_path)

        self.assertEqual(ds.sizes, arr2.sizes)
        np.testing.assert_array_equal(ds.coords["t"].values, arr2.coords["t"].values)
        np.testing.assert_array_equal(ds.data.values, arr2.values)
        
    def test_float_io(self):
        """Test reading/writing standard Float data."""
        # Create float32 data
        data = np.random.randn(10, 5).astype(np.float32)
        arr = xr.DataArray(
            data, 
            dims=("t", "x"),
            coords={"t": np.arange(10), "x": np.arange(5)}
        )

        # Write
        mx.xarray_to_rsf(arr, self.rsf_path)
        
        # Verify Header (m8r should see it as float)
        inp = m8r.Input(self.rsf_path)
        # RSF often reports this as "native_float" or similar
        self.assertIn("float", inp.string("data_format"))
        inp.close()

        # Read Back
        da_read = mx.rsf_to_xarray(self.rsf_path)
        
        # Assertions
        self.assertEqual(da_read.dtype, np.float32)
        np.testing.assert_allclose(da_read.values, arr.values)

    def test_int_io(self):
        """Test reading/writing Integer data."""
        # Create int32 data
        data = np.random.randint(-100, 100, size=(10, 5)).astype(np.int32)
        arr = xr.DataArray(
            data, 
            dims=("t", "x"),
            coords={"t": np.arange(10), "x": np.arange(5)}
        )

        # Write
        mx.xarray_to_rsf(arr, self.rsf_path)
        
        # Verify Header
        inp = m8r.Input(self.rsf_path)
        # RSF should report this as "native_int"
        self.assertIn("int", inp.string("data_format"))
        inp.close()

        # Read Back
        da_read = mx.rsf_to_xarray(self.rsf_path)
        
        # Assertions
        self.assertEqual(da_read.dtype, np.int32)
        np.testing.assert_array_equal(da_read.values, arr.values)

    def test_complex_io(self):
        """Test reading/writing Complex data."""
        # Create complex data
        nt, nx = 10, 5
        real_part = np.random.randn(nt, nx).astype(np.float32)
        imag_part = np.random.randn(nt, nx).astype(np.float32)
        data = (real_part + 1j * imag_part).astype(np.complex64)
        
        arr = xr.DataArray(
            data, 
            dims=("t", "x"),
            coords={"t": np.arange(nt), "x": np.arange(nx)}
        )

        # Write
        mx.xarray_to_rsf(arr, self.rsf_path)
        
        # Verify Header
        inp = m8r.Input(self.rsf_path)
        # RSF should report this as "native_complex"
        self.assertIn("complex", inp.string("data_format"))
        inp.close()

        # Read Back
        da_read = mx.rsf_to_xarray(self.rsf_path)
        
        # Assertions
        self.assertTrue(np.iscomplexobj(da_read))
        self.assertEqual(da_read.dtype, np.complex64)
        
        # Check values (real and imag parts)
        np.testing.assert_allclose(da_read.values, arr.values)  

    def test_dask_chunks(self):
        """Test Dask chunking functionality."""
        # Create large array
        dims = (100, 100, 100, 5)
        size = np.prod(dims)
        arr = xr.DataArray(
            np.arange(size, dtype=np.float32).reshape(dims),
            dims=("t","x","y","z"),
            coords={
                "t" : np.arange(100),
                "x" : np.arange(100),
                "y" : np.arange(100),
                "z" : np.arange(5)
            }
        )
        
        mx.xarray_to_rsf(arr, self.rsf_path)
        expected_mean = np.mean(np.arange(size, dtype=np.float32))
        del arr

        # Helper to time and check mean
        def check_chunks(chunk_config):
            start_time = time.time()
            da = mx.rsf_to_xarray(self.rsf_path, chunks=chunk_config)
            computed_mean = da.mean().compute().values
            duration = time.time() - start_time
            
            print(f"Chunks {chunk_config}: Mean={computed_mean}, Time={duration:.4f}s")
            self.assertTrue(np.isclose(computed_mean, expected_mean), 
                            f"Mean mismatch with chunks {chunk_config}")

        # Test various chunk configurations
        check_chunks((2, 3, 5, 1))
        # check_chunks((1, 1, 1, 1)) # Uncomment if you want to test no chunking
        check_chunks("auto")

    def test_sfwindow_xarray(self):
        """Test Madagascar sfwindow program on Xarray input."""
        # Create synthetic seismic data (Time, Offset)
        nt, nx = 100, 10
        dt, dx = 0.004, 25.0
        ot, ox = 0.0, 100.0
        
        data = np.zeros((nt, nx), dtype=np.float32)
        for ix in range(nx):
            shift = int(ix * 2) 
            if shift < nt:
                data[shift:, ix] = np.arange(nt - shift)

        arr = xr.DataArray(
            data, 
            dims=("t", "x"),
            coords={
                "t": np.arange(nt) * dt + ot,
                "x": np.arange(nx) * dx + ox
            }
        )

        # Subsampling (j1=2)
        res_sub = m8r.Filter('sfwindow')(j1=2).apply(arr)
        
        self.assertEqual(res_sub.sizes['t'], nt // 2)
        self.assertEqual(res_sub.sizes['x'], nx)
        new_dt = res_sub.coords['t'].values[1] - res_sub.coords['t'].values[0]
        self.assertTrue(np.isclose(new_dt, dt * 2))
        np.testing.assert_allclose(res_sub.values, arr.values[::2, :])

        # Windowing (f1, n1)
        f1, n1 = 10, 20
        res_win = m8r.Filter('sfwindow')(f1=f1, n1=n1).apply(arr)
        
        self.assertEqual(res_win.sizes['t'], n1)
        expected_o1 = ot + f1 * dt
        self.assertTrue(np.isclose(res_win.coords['t'].values[0], expected_o1))
        np.testing.assert_allclose(res_win.values, arr.values[f1:f1+n1, :])

        # Trace Selection (n2)
        n2 = 3
        res_trc = m8r.Filter('sfwindow')(n2=n2).apply(arr)
        
        self.assertEqual(res_trc.sizes['x'], n2)
        np.testing.assert_array_equal(res_trc.values, arr.values[:, :n2])

        # Piping (j1=2 | n1=10)
        op = m8r.Filter('sfwindow')(j1=2) | m8r.Filter('sfwindow')(n1=10)
        res_pipe = op.apply(arr)
        
        self.assertEqual(res_pipe.sizes['t'], 10)
        new_dt_pipe = res_pipe.coords['t'].values[1] - res_pipe.coords['t'].values[0]
        self.assertTrue(np.isclose(new_dt_pipe, dt * 2))

    def test_sfwindow_numpy(self):
        """Test Madagascar sfwindow program on Numpy input (Legacy Support)."""
        nt, nx = 100, 10
        
        # Create Data in (Offset, Time) format for Numpy
        data = np.zeros((nx, nt), dtype=np.float32)
        for ix in range(nx):
            shift = int(ix * 2) 
            if shift < nt:
                data[ix, shift:] = np.arange(nt - shift)

        arr = data

        # Subsampling (j1=2) 
        res_sub = m8r.Filter('sfwindow')(j1=2).apply(arr)
        self.assertEqual(res_sub.shape[1], nt // 2)
        self.assertEqual(res_sub.shape[0], nx)
        np.testing.assert_allclose(res_sub, arr[:, ::2])

        # Windowing (f1, n1) 
        f1, n1_win = 10, 20
        res_win = m8r.Filter('sfwindow')(f1=f1, n1=n1_win).apply(arr)
        self.assertEqual(res_win.shape[1], n1_win)
        np.testing.assert_allclose(res_win, arr[:, f1:f1+n1_win])

        # Trace Selection (n2) 
        n2_cut = 3
        res_trc = m8r.Filter('sfwindow')(n2=n2_cut).apply(arr)
        self.assertEqual(res_trc.shape[0], n2_cut)
        np.testing.assert_array_equal(res_trc, arr[:n2_cut, :])

        # Piping
        op = m8r.Filter('sfwindow')(j1=2) | m8r.Filter('sfwindow')(n1=10)
        res_pipe = op.apply(arr)
        self.assertEqual(res_pipe.shape[1], 10)
        
        expected_pipe = arr[:, ::2][:, :10]
        np.testing.assert_allclose(res_pipe, expected_pipe)
        
    def test_auxiliary_file(self):
        """Test passing an xarray as a keyword argument (e.g. velocity=Vp)."""
        # Create a dummy Velocity model (Vp)
        vp_arr = xr.DataArray(
            np.ones((10, 10), dtype=np.float32) * 2.0, 
            dims=("z", "x"),
            coords={"z": np.arange(10), "x": np.arange(10)}
        )
        
        # Create Data to process
        data_arr = xr.DataArray(
            np.zeros((10, 10), dtype=np.float32), 
            dims=("z", "x"),
            coords={"z": np.arange(10), "x": np.arange(10)}
        )

        # Run a command that takes an aux file. 
        math_op = m8r.Filter('sfmath')(output="input+vel", vel=vp_arr)
        
        res = math_op.apply(data_arr)
        
        # Result should be 0 + 2.0 = 2.0
        self.assertEqual(res.mean(), 2.0)

if __name__ == "__main__":
    unittest.main(verbosity=1)