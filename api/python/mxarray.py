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

import os
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

    rsf_type = getattr(f, 'type', None)
    if not rsf_type:
        fmt = f.string("data_format")
        if fmt and "complex" in fmt: rsf_type = 'complex'
        elif fmt and "int" in fmt: rsf_type = 'int'
        else: rsf_type = 'float'
    
    dtype_map = {
        'float':   np.float32,
        'int':     np.int32,
        'complex': np.complex64,
        'uchar':   np.uint8,
        'char':    np.int8
    }
    dtype = dtype_map.get(rsf_type, np.float32)
    # -------------------------------
    
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
    
    if np.issubdtype(ds.dtype, np.complexfloating):
        rsf_type = 'complex'
        out_dtype = np.complex64
        fmt_str = "native_complex"
    elif np.issubdtype(ds.dtype, np.integer):
        rsf_type = 'int'
        out_dtype = np.int32
        fmt_str = "native_int"
    else:
        rsf_type = 'float' 
        out_dtype = np.float32
        fmt_str = "native_float"
    # -------------------------------

    data = ds.values
    data = data.transpose(*reversed(range(data.ndim)))
    
    dims = ds.dims
    out = m8r.Output(outpath)

    # Set Type
    out.settype(rsf_type)
    out.put("data_format", fmt_str)
    
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
        
    out.write(data.astype(out_dtype))
    out.close()
    
def rsf_to_xarrayds(path, chunks = "auto"):
    """Convert an RSF file to an xarray Dataset."""
    
    da = rsf_to_xarray(path, chunks=chunks)
    ds = da.to_dataset(name="data")
    
    return ds

## Monkey patching m8r.Filter to handle xarray inputs/outputs

def _patched_setcommand(self, kw, args=[]):
    """
    Patched version of Filter.setcommand to handle auxiliary xarrays.
    Example: Filter('sfprog')(velocity=my_xarray)
    """
    if not hasattr(self, '_mx_aux_refs'):
        self._mx_aux_refs = []

    for key, val in list(kw.items()):
        if isinstance(val, xr.DataArray):
            tmp_name = m8r.Temp()
            xarray_to_rsf(val, tmp_name)

            f_obj = m8r.File(tmp_name, temp=True)
            self._mx_aux_refs.append(f_obj)
            
            kw[key] = str(tmp_name)

    new_args = []
    for val in args:
        if isinstance(val, xr.DataArray):
            tmp_name = m8r.Temp()
            xarray_to_rsf(val, tmp_name)
            f_obj = m8r.File(tmp_name, temp=True)
            self._mx_aux_refs.append(f_obj)
            new_args.append(str(tmp_name))
        else:
            new_args.append(val)
            
    # Call the original function with strings/files only
    original_func = getattr(self.__class__, '_original_setcommand_func', None)
    if not original_func:
        # Fallback if class attribute missing
        real_module = m8r.wrapped if hasattr(m8r, 'wrapped') else m8r
        original_func = real_module.Filter._original_setcommand_func
        
    return original_func(self, kw, new_args)


def _patched_apply(self, *srcs):
    """
    Patched version of Filter.apply to handle xarray inputs/outputs.
    """
    
    # Check if any input is an xarray
    if not any(isinstance(s, xr.DataArray) for s in srcs):
        # Fallback to original
        if hasattr(self.__class__, '_original_apply_func'):
            return self.__class__._original_apply_func(self, *srcs)
        real_module = m8r.wrapped if hasattr(m8r, 'wrapped') else m8r
        return real_module.Filter._original_apply_func(self, *srcs)

    # Handle xarray Input
    clean = []
    rsf_inputs = []
    try:
        for s in srcs:
            if isinstance(s, xr.DataArray):
                tmp = m8r.Temp()
                xarray_to_rsf(s, tmp)
                rsf_inputs.append(str(tmp))
                clean.extend([str(tmp), str(tmp)+'@'])
            else:
                rsf_inputs.append(str(s))

        out_file = m8r.Temp()
        
        # Split command for correct pipe handling
        first, pipe_char, rest = self.command.partition('|')
        cmd_parts = [first]
        
        if len(rsf_inputs) > 0:
            cmd_parts.append(f"< {rsf_inputs[0]}")
            if len(rsf_inputs) > 1:
                cmd_parts.extend(rsf_inputs[1:])
        
        if pipe_char:
            cmd_parts.append(pipe_char)
            cmd_parts.append(rest)
            
        cmd_parts.append(f"> {out_file}")
        
        full_cmd = " ".join(cmd_parts)
        
        if os.system(full_cmd) != 0:
            raise RuntimeError(f"Command failed: {full_cmd}")

        if self.plot: 
            return m8r.Vplot(out_file, temp=True)
            
        res = rsf_to_xarray(out_file)
        clean.extend([str(out_file), str(out_file)+'@'])
        return res

    finally:
        for f in clean:
            if os.path.exists(f): 
                try: os.unlink(f)
                except: pass
                
if hasattr(m8r, 'wrapped'):
    real_m8r_module = m8r.wrapped
else:
    real_m8r_module = m8r

if not hasattr(real_m8r_module.Filter, '_original_apply_func'):
    real_m8r_module.Filter._original_apply_func = real_m8r_module.Filter.apply
    real_m8r_module.Filter.apply = _patched_apply

if not hasattr(real_m8r_module.Filter, '_original_setcommand_func'):
    real_m8r_module.Filter._original_setcommand_func = real_m8r_module.Filter.setcommand
    real_m8r_module.Filter.setcommand = _patched_setcommand