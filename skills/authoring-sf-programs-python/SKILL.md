---
name: authoring-sf-programs-python
description: Use when authoring a Madagascar sf* main program in Python (fast prototyping with numpy).
---

## When to use

Load this skill when writing a new `sf<name>` program in Python. Python is the fastest iteration path: no compilation step, numpy-native array operations, and easy debugging. New Python programs go in `user/<youruser>/M<name>.py`.

This skill covers the Python-specific layer only. Always also load:

- `../authoring-sf-programs/SKILL.md` — file naming conventions, self-documentation requirements, parameter style, and build integration shared by every language.

Use this skill instead of the C skill whenever:
- You are prototyping an algorithm and want fast iteration.
- Your operation is naturally expressed as numpy array math.
- You do not need raw C performance (inner loops in pure Python can be 10–100x slower than C; use numpy vectorization to avoid that).

---

## Skeleton

This skeleton has been smoke-tested against a real `sfspike`-generated RSF file and runs correctly. The key difference from the plan's template: the correct import is `rsf.api as rsf` (the installed module name), not `import m8r`, and `rsf.Output()` takes only a tag argument — shape is inherited automatically from the first input opened.

```python
#!/usr/bin/env python3
'''One-line description of what this program does.'''

import sys
import numpy as np
import rsf.api as rsf

par = rsf.Par()
fin  = rsf.Input()        # stdin by default
fout = rsf.Output()       # stdout by default; inherits format from fin

n1 = fin.int("n1")        # fastest axis length (from RSF header)
n2 = fin.size(1)          # product of all axes beyond n1 (number of traces)

factor = par.float("factor", 1.0)  # scale factor [default 1.0]

trace = np.zeros(n1, dtype='f')

for i2 in range(n2):
    fin.read(trace)
    trace *= factor
    fout.write(trace)

sys.exit(0)
```

### Bulk-read variant (whole array at once)

When the algorithm needs the full dataset in memory, read everything at once:

```python
#!/usr/bin/env python3
'''One-line description.'''

import sys
import numpy as np
import rsf.api as rsf

par = rsf.Par()
fin  = rsf.Input()
fout = rsf.Output()

n1 = fin.int("n1")
n2 = fin.int("n2") or 1

factor = par.float("factor", 1.0)  # scale factor [default 1.0]

data = np.zeros((n2, n1), dtype=np.float32)
fin.read(data)

data *= factor

fout.write(data)
sys.exit(0)
```

### Correction from the plan's template

The plan template showed `m8r.Output("out", fin)`. This is wrong in two ways:

1. Prefer `import rsf.api as rsf` over `import m8r`. Both work (the installed `m8r.py` is bit-identical to `rsf/api.py`), but the self-doc scraper's I/O-recognition regex (`inpout['python']` in `framework/rsf/doc.py`) is anchored to `rsf.Input` / `rsf.Output` calls. Using `import m8r` would lose the auto-generated `< in.rsf > out.rsf` synopsis in `sfdoc` output. Note that `api/python/test/clip.py` and `api/python/test/afdm.py` both use `import m8r` — they still work, but their sfdoc entries are less informative.
2. `rsf.Output(tag='out', data_format=None)` accepts only a tag and an optional format string — there is no second positional argument for a source file. Shape inheritance is automatic: when writing begins, the internal `_RSF.fileflush()` copies header parameters from the first opened input.

---

## Self-doc format

The scraper in `framework/rsf/doc.py` uses this regex for Python files (`comment['python']`):

```python
comment['python'] = re.compile(
    r'[^\'\"]*([\'\"]+)(?P<comment>[^\'\"].+?)\1',
    re.DOTALL
)
```

This matches the **first string literal** in the file — specifically the first `'...'` or `"""..."""` that appears before any other quoted string. In practice: put a module-level docstring immediately after the shebang line:

```python
#!/usr/bin/env python3
'''Scale input data by a scalar factor.'''
```

or a triple-quoted version:

```python
#!/usr/bin/env python3
"""
Scale input data by a scalar factor.

Longer explanation on subsequent lines.
"""
```

The first line of the docstring becomes the program's one-line description shown in `sfdoc sfname`. Subsequent lines become the comments section.

Real example from `api/python/test/afdm.py`: that file has no leading docstring, so its `sfdoc` description is empty. Real example from `user/sbader/Menergy.py`:

```python
#!/usr/bin/env python3
'''
Estimate energy of input

E(t) = \sum\limits_{k=(t-\frac{R}{2})}^{(t+\frac{R}{2})}A(k)^2
'''
```

The parameter scraper regex for Python (`param['python']`) captures `par.bool(...)`, `par.int(...)`, `par.float(...)`, and `par.string(...)` calls along with optional trailing `# [range] description` comments:

```python
param['python'] = re.compile(
    r'par\.(?P<type>bool|int|float|string)'
    r'\s*\(\s*[\"\'](?P<name>\w+)[\"\']\s*'
    r'(?:\,\s*(?P<default>[^\)]+))?\)'
    r'(?:\s*\#\s*(?P<range>[\[][^\]]+[\]])?\s*'
    r'(?P<desc>[^#\n]+\S))?'
)
```

So a parameter line like:

```python
factor = par.float("factor", 1.0)  # [0,inf] scale factor applied to every sample
```

scrapes as: type=float, name=factor, default=1.0, range=[0,inf], desc="scale factor applied to every sample".

---

## m8r.Input / m8r.Output

Both classes live in `api/python/m8r.py` (installed as `rsf/api.py`). They extend `_File`, which extends `File`.

### `rsf.Input(tag='in')`

Opens an RSF file for reading.

- `tag='in'` means standard input (the default, used when the program reads from a pipe).
- Any other tag string (e.g., `'vel'`, `'ref'`) opens the file named by the command-line argument `vel=some.rsf`.

From `afdm.py`:
```python
Fw = rsf.Input()         # stdin — the wavelet file
Fv = rsf.Input("vel")    # opened via vel= on the command line
Fr = rsf.Input("ref")    # opened via ref= on the command line
```

**Header reads** — extract axis metadata from the RSF header:

```python
n1 = fin.int("n1")           # returns int or None
d1 = fin.float("d1")         # returns float or None
n1 = fin.int("n1", default=100)   # returns default if key absent
label = fin.string("label1") # returns str or None

# Axis convenience method (returns dict with n, d, o, l, u):
ax = fin.axis(1)
nt = ax['n']; dt = ax['d']; ot = ax['o']
```

**Reading data:**

```python
trace = np.zeros(n1, dtype='f')
fin.read(trace)       # fills trace in-place (trace-by-trace loop pattern)

n2 = fin.size(1)      # product of all axes from axis 2 onward
data = np.zeros((n2, n1), dtype='f')
fin.read(data)        # bulk read; shape (n2, n1)
```

`fin.close()` is optional — the destructor closes automatically.

### `rsf.Output(tag='out', data_format=None)`

Opens an RSF file for writing.

- `tag='out'` means standard output (the default).
- Shape and format are inherited automatically from the first `rsf.Input()` opened in the same program; you do not pass a source file.
- Use `.put()` to override or add header values before any `.write()` call.

**Writing metadata:**

```python
fout.put("n1", n1)
fout.put("d1", 0.004)
fout.put("label1", "Time")
fout.put("unit1", "s")

# Axis convenience:
fout.putaxis(ax, 1)   # copies n, d, o, label, unit from axis dict
```

**Writing data:**

```python
fout.write(trace)       # write a 1D numpy array
fout.write(data)        # write a 2D (or ND) numpy array — flattened internally
```

The `write()` method calls `np.reshape(data.astype(np.float32), (data.size,))` before writing, so the shape of the array passed does not need to match the file shape exactly — the total element count must match.

**Closing:**

```python
fout.close()   # optional; destructor closes automatically
```

---

## m8r.Par

`rsf.Par(argv=sys.argv)` parses command-line arguments of the form `key=value`. A default instance is created automatically on module import; call `rsf.Par()` explicitly in your program to initialize it properly.

### Scalar parameters

```python
par = rsf.Par()

i   = par.int("niter")           # returns int or None
i   = par.int("niter", 10)       # returns 10 if niter= not on command line
f   = par.float("eps", 0.01)     # float with default
b   = par.bool("verb", False)    # bool: True if verb=y/Y/1; False if n/N/0
s   = par.string("mode", "fwd")  # string (strips surrounding quotes)
```

Booleans on the command line: `verb=y` or `verb=1` → True; `verb=n` or `verb=0` → False.

### List (array) parameters

```python
ints   = par.ints("k1", 3)            # read 3 ints: k1=1,5,9
floats = par.floats("scale", 2)       # read 2 floats: scale=1.0,2.0
bools  = par.bools("flags", 4)        # read 4 bools
```

If fewer values are given than requested, the last value is repeated to fill.

### Program name

```python
name = par.getprog()   # returns sys.argv[0]
```

### Parameter file

Pass `par=somefile.par` on the command line to read `key=value` pairs from a file; they merge with command-line arguments.

---

## Numpy shape conventions

Madagascar numbers axes starting from 1, with **axis 1 being the fastest-varying** (stored contiguously in memory, analogous to C's inner loop). Numpy uses C order (row-major), where the **last index is fastest-varying**.

Therefore, when you allocate a numpy array to hold Madagascar data, the axis ordering is **reversed**:

| Madagascar | Numpy shape |
|------------|-------------|
| n1=100 samples/trace | last dimension: `(..., 100)` |
| n2=50 traces | second-to-last: `(50, 100)` |
| n3=10 shots | `(10, 50, 100)` |

**Example**: a 3D dataset with n1=500, n2=200, n3=10:

```python
n1 = fin.int("n1")   # 500 — fastest (time samples)
n2 = fin.int("n2")   # 200 — traces
n3 = fin.int("n3")   # 10  — shots
data = np.zeros((n3, n2, n1), dtype=np.float32)
fin.read(data)
# data[ishot, itrace, isample] — numpy index order
```

### Common transpose patterns

When an algorithm expects "time along rows" (C convention: time last), your layout is already correct. When an external library expects "columns are samples" (Fortran convention), transpose:

```python
# Madagascar stores (n2, n1) — n2 rows, n1 samples per row
data = np.zeros((n2, n1), dtype='f')
fin.read(data)

# For a library expecting (n1, n2) — samples along rows:
data_T = data.T   # view, no copy

# Or use numpy's built-in where a copy is needed:
data_fortran = np.asfortranarray(data)
```

When writing back, remember to restore the Madagascar axis order:

```python
result = some_lib_call(data.T)   # result shape is (n1, n2)
fout.write(result.T)             # write (n2, n1) — Madagascar order
```

### `fin.shape()` and axis reversal

The `File.shape()` method returns a tuple already reversed for numpy:

```python
sh = fin.shape()   # e.g. (10, 50, 100) for n1=100, n2=50, n3=10
data = np.zeros(sh, dtype='f')
fin.read(data)
```

---

## SConstruct integration

### Declaring Python programs

In `user/<youruser>/SConstruct`, list each Python main program name (without the `M` prefix and without `.py`) in the `pyprogs` variable:

```python
import os, sys, re
sys.path.append('../../framework')
import bldutil

pyprogs = 'myscale mystack'   # Mscale.py and Mstack.py

try:
    Import('env root pkgdir bindir libdir incdir')
    env = env.Clone()
except:
    env = bldutil.Debug()
    root = None

if root: # no compilation, just rename
    pymains = Split(pyprogs)
    exe = env.get('PROGSUFFIX','')
    for prog in pymains:
        binary = os.path.join(bindir,'sf'+prog+exe)
        env.InstallAs(binary,'M'+prog+'.py')
        env.AddPostAction(binary,Chmod(str(binary),0o755))
```

### What the build does

When you run `scons` in the top-level source tree, the build system processes `user/<youruser>/SConstruct`. For each name in `pyprogs`, the inline install loop does:

1. Copies `M<name>.py` directly to `$RSFROOT/bin/sf<name>` (no compilation) via `env.InstallAs`.
2. Sets the file executable (`chmod 0o755`) via `env.AddPostAction`.

There is **no shell wrapper** generated for Python programs. The Python file itself becomes the installed binary, called directly by the system Python interpreter via the shebang line (`#!/usr/bin/env python3`). This is different from C programs (which are compiled ELF binaries) but identical in usage from the user's perspective.

### Using Python programs in a Flow

Once installed, use `sf<name>` exactly like any other Madagascar program:

```python
Flow('out', 'in', 'sfmyscale factor=2.0')
Flow('out', ['in', 'vel'], 'sfmystack vel=${SOURCES[1]}')
```

For testing in `user/<youruser>/`, you can run the Python file directly without installing:

```bash
python3 Mmyscale.py factor=2.0 < in.rsf > out.rsf
```

---

## Worked references

- **`api/python/test/clip.py`** — minimal example. Uses a trace-by-trace loop with `fin.size(1)` to count traces; calls `par.float("clip")` (required, no default) and clips with `numpy.clip`. ~25 lines.

- **`api/python/test/afdm.py`** — richer example: acoustic finite-difference modeling. Opens three inputs (`Fw`, `Fv`, `Fr`) and one output (`Fo`); uses `fin.axis()` / `fout.putaxis()` for full axis metadata; reads full arrays then loops in time; writes one time snapshot per iteration. ~70 lines.

- **`user/sbader/Menergy.py`** — real user program. Triple-quoted docstring for self-doc; uses `rsf.api as rsf`; reads a 1D dataset and computes a rolling energy estimate. Shows that `rsf.Output()` can be opened after `rsf.Input()` is already open (the output inherits the input's format automatically).

- **`user/sbader/Mreplace.py`** — another simple 1D user program. Good minimal reference for the `rsf.api` import pattern.

- **`user/godwinj/Mpysvd.py`** — real user program that wraps SciPy SVD. Demonstrates the recommended try/except guard around the `rsf.api` import: if the Python API or NumPy/SciPy are absent the program prints a clear error rather than a raw `ImportError`. Use this pattern whenever your program has optional heavyweight dependencies.

---

## Building and testing

### Local run (no install needed)

```bash
cd user/<youruser>/
python3 Mmyscale.py factor=2.0 < in.rsf > out.rsf
sfattr < out.rsf
```

### Install via scons

```bash
cd /path/to/RSFSRC
scons user/<youruser>/
```

After a successful build, `$RSFROOT/bin/sfmyscale` is the Python file itself (executable). Test it:

```bash
sfmyscale factor=2.0 < in.rsf > out.rsf
sfattr < out.rsf
```

### Self-doc check

```bash
sfdoc sfmyscale
```

If the program description is empty, the leading docstring was not found. Verify it is a bare string literal on the second line (right after the shebang), not a comment.

### Regression test

Add a test flow in the SConstruct:

```python
Flow('test_out', 'test_in', 'sfmyscale factor=2.0')
```

Then `scons test_out` runs the program as part of the build. For automated comparison, use `sfmath`:

```python
Flow('diff', ['test_out', 'reference'], 'sfmath x=${SOURCES[0]} y=${SOURCES[1]} output="x-y" | sfattr want=max')
```

### Common pitfalls

- `ImportError: No module named rsf.api` — source the Madagascar env: `source $RSFROOT/etc/env.sh`.
- Empty `sfdoc` description — move the docstring to line 2 right after the shebang; it must be a string literal, not a `#` comment.
- `TypeError: Unsupported type in put` — numpy scalars (int64, float64) may need explicit casting: `fout.put("n1", int(n1))`.
- Wrong output dimensions — verify `data.size == n1 * n2 * ...` before calling `fout.write(data)`.
