import numpy as np
import warnings, re, os, io, datetime, socket
from .utils import _check_input_source, _str_match_re, _get_datapath
from .version import __version__

RSFHSPLITER = b"\x0c\x0c\x04"
io.BytesIO()
def read_rsf(file, order='F'):
    """
    Read RSF file and return (data, header) or None.

    Parameters
    ----------
    file : str or file-like object
        The RSF file to read.

    Returns
    -------
    list [ndarray, dict, str]
        A list containing the read data, header, and history information.
    order : str
        The order in which to read the data (default is 'F' for Fortran-style).
    """
    try:
        close_after = isinstance(file, str)

        file_fp = _check_input_source(file, 'rb')
        if file_fp is None:
            return None

        header = {}
        data = None

        # Read header
        buf = bytearray()
        while True:
            chunk = file_fp.read(1)
            if not chunk:
                break
            buf.extend(chunk)
            if buf.endswith(RSFHSPLITER):
                break

        header_text = buf.rstrip(RSFHSPLITER).decode("utf-8", errors="ignore")
        
        header = _str_match_re(header_text)

        # Format conversion
        for k, v in list(header.items()):
            if re.fullmatch(r"n[1-9]", k) or k == "esize":
                try:
                    header[k] = int(v)
                except ValueError:
                    pass
            elif re.fullmatch(r"[od][1-9]", k):
                try:
                    header[k] = float(v)
                except ValueError:
                    pass

        # data source
        in_val = header.get("in", None)
        if in_val is None:
            warnings.warn("'in' key not found in RSF header")
            return None

        if in_val == "stdin":
            data_file = file_fp
        else:
            data_file = _check_input_source(in_val, 'rb')
            if data_file is None:
                warnings.warn(f"Data file not accessible: {in_val}")
                return None

        # shape
        shape = []
        for i in range(1, 10):
            key = f"n{i}"
            if key in header:
                shape.append(int(header[key]))
            else:
                break
        if not shape:
            warnings.warn("No n# keys found for shape")
            return None

        # data_format
        fmt = header.get("data_format", None)
        if fmt is None:
            warnings.warn("'data_format' key not found")
            return None
        try:
            fmt_A, fmt_B = fmt.split("_", 1)
        except ValueError:
            warnings.warn(f"Invalid data_format: {fmt}")
            return None

        if fmt_A not in ("native", "ascii", "xdr"):
            warnings.warn(f"Unsupported format type: {fmt_A}")
            return None
        if fmt_B not in ("int", "float", "complex"):
            warnings.warn(f"Unsupported data type: {fmt_B}")
            return None

        dtype_map = {
            "int": np.int32,
            "float": np.float32,
            "complex": np.complex64
        }
        dtype = dtype_map[fmt_B]

        # 读取数据
        if fmt_A == "ascii":
            ascii_text = data_file.read().decode("utf-8", errors="ignore").strip()
            parts = ascii_text.split()
            if fmt_B == "complex":
                def parse_complex(s):
                    s = s.replace("i", "j")
                    return complex(s)
                arr = np.array([parse_complex(p) for p in parts], dtype=dtype)
            else:
                arr = np.array([float(p) for p in parts], dtype=dtype)
        else:
            arr = np.frombuffer(data_file.read(), dtype=dtype)
            if fmt_A == "xdr":
                swapped = arr.byteswap()
                if hasattr(swapped, 'newbyteorder'):
                    arr = swapped.newbyteorder()
                else: arr = arr.view(swapped.dtype.newbyteorder())

        arr = arr.reshape(shape, order='F')

        if in_val != "stdin":
            data_file.close()
        if close_after:
            file_fp.close()

        return [arr, header, header_text]

    except Exception as e:
        warnings.warn(f"Error reading RSF: {e}")
        return None



def write_rsf(arr: np.ndarray, file, header={}, history='', out=None, form="native", fmt="%f"):
    """
    Write RSF file with given header and data.

    Parameters:
    arr : ndarray
        The data array to write.
    file : str or file-like object
        The output file (header) or file-like object.
    header : dict
        The header information to write.
    out : str, optional
        Data bytes output file or file-like object.
    form : str, optional
        The data format (default is "native").
        native: little-endian
        xdr: network (big-endian) byte order
        ascii: plain text
    """
    close_after = isinstance(file, str)

    outheader = {}
    outheader.update(header if isinstance(header, dict) else {})
    if not isinstance(arr, np.ndarray):
        raise TypeError(f"Expected ndarray, got {type(arr)}")

    file_fp = _check_input_source(file, 'wb')

    if file_fp is None:
        raise ValueError(f"Cannot open file: {file}")
    
    if out is None:
        if isinstance(file, str):
            fname = os.path.basename(file)
            fname = os.path.join(_get_datapath(), fname + '@')
            out_fp = _check_input_source(fname, 'wb')
            outheader["in"] = os.path.abspath(fname)
        else:
            out_fp = file_fp
        
    elif out == 'stdout':
        out_fp = file_fp
    elif isinstance(out, str) or isinstance(out, io.IOBase):
        out_fp = _check_input_source(out, 'wb')
        if out_fp is None:
            raise ValueError(f"Cannot open output file: {out}")
        if isinstance(out, str):
            outheader["in"] = os.path.abspath(out)
        else:
            outheader["in"] = out

    dtype = arr.dtype.name
    dtype = ''.join([c for c in dtype if c.isalpha()])
    outheader.update({"data_format": f"{form}_{dtype}"})
    # cast outheader to string
    header_str = ""
    for key, value in outheader.items():
        if isinstance(value, str):
            if key not in ["in","data_format","esize","out"]: value = f'"{value}"'
            header_str += f"{key}={value}\n"
        else:
            header_str += f"{key}={str(value)}\n"

    # Write header
    try:
        uname = os.getlogin()
        hostname = socket.gethostname()
    except:
        uname = "unknown"
        hostname = "unknown"
    banner = f"RSFPY_{__version__}\t{os.getcwd()}\t{uname}@{hostname}\t{datetime.datetime.now().strftime('%a %b %d %H:%M:%S %Y')}"
    file_fp.write(history.encode('utf-8') + b"\n\n")
    file_fp.write(banner.encode('utf-8') + b"\n")
    file_fp.write(header_str.encode('utf-8') + b"\n\n")
    if out is None or out == 'stdout': file_fp.write(RSFHSPLITER)

    # Write data
    if form == "ascii":
        np.savetxt(out_fp, arr, fmt=fmt)
    elif form == "native":
        out_fp.write(arr.tobytes())
    elif form == "xdr":
        out_fp.write(arr.byteswap().tobytes())

    if close_after:
        file_fp.close()
