---
name: authoring-sf-programs-matlab
description: Use when authoring a Madagascar sf* main program in MATLAB.
---

## When to use

Use this skill when a researcher already has MATLAB installed (licensed) and wants
to write an `sf<name>` program as a `.m` script, leveraging existing MATLAB toolboxes
(e.g. Signal Processing Toolbox, Optimization Toolbox) that would be cumbersome to
replicate in C or Python.

**File naming**: the MATLAB main script must be named `M<name>.m` — same `M` prefix
convention used by all Madagascar languages. There is no `sf` prefix in the source file;
the installed binary is called `sf<name>`.

**Important — compile-time requirement**: the MEX backing library (`rsf_create.mex*`,
`rsf_read.mex*`, etc.) is compiled at `scons` configure time only when MATLAB is
detected on the build host. See `api/matlab/SConstruct`: the entire MEX build is gated
on `'matlab' in env.get('API',[])`. Without a licensed MATLAB on the build machine the
API produces nothing and MATLAB programs cannot run.

**Octave overlap**: the `.m` extension is shared with Octave. The build system selects
Matlab vs. Octave at configure time via the `API` variable; they are mutually exclusive
at build time. For Octave programs see `authoring-sf-programs-octave/SKILL.md`.

Always load the shared skill alongside this one:
[`skills/authoring-sf-programs/`](../authoring-sf-programs/SKILL.md) — file naming,
self-documentation, parameter conventions, error handling, and testing conventions that
apply to every sf-program regardless of language.

---

## Skeleton

The canonical pattern is taken from `api/matlab/test/clip.m`. A MATLAB sf-program
is a plain function (not a script) that receives RSF filenames as arguments:

```matlab
function M<name>(in, out, param1)
%<NAME>  One-sentence synopsis (appears in sfdoc output).
%
% par1=default  description of parameter

% 1. Query dimensions of the input file.
dims = rsf_dim(in);          % returns column vector [n1; n2; ...]
n1   = dims(1);
n2   = prod(dims(2:end));

% 2. Read a named parameter from the input header.
clip = rsf_par(in, 'clip', 'f', 1.0);  % type 'f'=float, 'i'=int, 'b'=bool

% 3. Create output header, inheriting geometry from input.
rsf_create(out, in);         % copy header: rsf_create(outfile, infile)
% -- OR, specify sizes explicitly:
% rsf_create(out, [n1; n2]);

% 4. Allocate a trace-length working buffer.
trace = zeros(1, n1);

% 5. Loop over traces, reading and writing sequentially.
for i2 = 1:n2
    rsf_read(trace, in, 'same');          % read next n1 samples
    % ... process trace ...
    trace(trace >  clip) =  clip;
    trace(trace < -clip) = -clip;
    rsf_write(trace, out, 'same');        % write next n1 samples
end
```

The `'same'` flag on `rsf_read` / `rsf_write` advances an internal file pointer so
successive calls walk through the binary data sequentially — essential for loop-based
trace processing.

For single-call whole-file I/O, prefer `rsf_read_all` / `rsf_write_all` (see cheat
sheet below).

---

## API cheat sheet

All functions are MEX-compiled C entry points. `nrhs` = number of right-hand side
(input) arguments; `nlhs` = number of left-hand side (output) arguments.

| Function | MATLAB call | nrhs | nlhs | Notes |
|---|---|---|---|---|
| `rsf_create` | `rsf_create(outfile, infile)` | 2 | 0 | Create RSF header on disk; copies geometry from `infile` (string) |
| `rsf_create` | `rsf_create(outfile, [n1;n2;...])` | 2 | 0 | Create header; sizes from column vector of doubles |
| `rsf_read` | `rsf_read(buf, infile)` | 2 | 0 | Read `numel(buf)` samples from `infile` into pre-allocated double array `buf` |
| `rsf_read` | `rsf_read(buf, infile, 'same')` | 3 | 0 | Same, but advance internal file pointer (sequential trace loop) |
| `rsf_read_all` | `[data,sz,d,o,lbl,unit] = rsf_read_all(file)` | 1 | 1–6 | Read entire file; returns data array plus optional size/delta/origin/label/unit vectors |
| `rsf_read_header` | `[sz,d,o,lbl,unit] = rsf_read_header(file)` | 1 | 1–5 | Read header only (no data); returns size vector plus optional delta/origin/label/unit |
| `rsf_write` | `rsf_write(buf, outfile)` | 2 | 0 | Write double array `buf` to `outfile` (sets n1/n2/... from `buf` dimensions) |
| `rsf_write` | `rsf_write(buf, outfile, 'same')` | 3 | 0 | Append to existing binary; does not rewrite header |
| `rsf_write_all` | `rsf_write_all(file, cmdargs, data)` | 3–7 | 0 | Write full RSF file; `cmdargs` is cell array (e.g. `{'out=stdout'}`); optional delta/origin/label/unit row vectors |
| `rsf_par` | `v = rsf_par(infile, 'name', 'f', default)` | 4 | 1 | Read scalar header parameter; type string: `'f'`=float, `'i'`/`'d'`=int, `'l'`/`'b'`=bool; returns `default` if absent |
| `rsf_dim` | `dims = rsf_dim(infile)` | 1 | 1 | Return column vector of axis lengths `[n1; n2; ...]` for the RSF file |

**Type details for `rsf_par`**: type string first character selects the branch —
`'f'` or `'r'` → float, `'i'` or `'d'` → integer, `'l'` or `'b'` → logical.
The return value is always a double scalar.

**Complex data**: `rsf_read` / `rsf_write` handle `SF_COMPLEX` files automatically
when `buf` / `data` is a MATLAB complex array (`mxIsComplex`).

**`m8r` MEX** (`m8r.c` → `m8r.mex*`): exposes `rsf_filter(cmd, data, deltas, origins,
labels, units)` — runs an arbitrary Madagascar binary as a filter on MATLAB data via
temp files. Useful for calling existing `sf*` programs from within a MATLAB session.

---

## Build integration

The MEX binaries are built by `api/matlab/SConstruct`, which is invoked automatically
by the top-level Madagascar `scons` when MATLAB is present:

```python
# api/matlab/SConstruct (excerpt)
if 'matlab' in env.get('API',[]):
    mex = env.get('MEX')          # path to mex compiler, detected at configure time
    suffix = env.get('MEXSUFFIX') # platform suffix, e.g. .mexa64 / .mexmaci64
    if mex:
        mexcom = mex + " CC=$CC CFLAGS='$CFLAGS $_CCCOMCOM -fPIC' ..."
        for inp in Split('m8r par dim read_header read read_all write write_all create'):
            ...
            mfile = env.Command(cfile+suffix, cfile+'.c', mexcom)
            if root:
                install = env.Install(libdir, mfile)
```

MATLAB must be detectable by the configure step for this block to execute. If MATLAB
is not installed, `env.get('API',[])` will not contain `'matlab'` and the entire block
is skipped — no MEX files are produced and MATLAB programs will fail with "undefined
function" errors at runtime.

**No `UserSconsTargets` for MATLAB**: the `bldutil.UserSconsTargets()` helper in
`framework/bldutil.py` supports `.c`, `.f90`, `.py`, and `.jl` attributes — not `.m`.
MATLAB programs are not automatically discovered or compiled by the standard user
directory `SConstruct` machinery.

**Deploying a new MATLAB program**: place `M<name>.m` somewhere on MATLAB's path (or
in the current directory) and ensure that the MEX binaries installed to `libdir` are
also on MATLAB's path. The simplest approach:

```matlab
addpath('/path/to/madagascar/lib')   % directory containing rsf_*.mex* files
addpath('/path/to/your/Mname.m')
M<name>(input_rsf, output_rsf, ...)
```

There is no automated install target for `.m` main programs analogous to the C or
Python install targets.

---

## Pointers

Files in `api/matlab/` with one-line descriptions:

| File | Description |
|---|---|
| `SConstruct` | SCons build script; compiles all `rsf_*.c` and `m8r.c` to MEX when `matlab` is in `API` |
| `rsf_create.c` | MEX: create an RSF output header from a filename+infile or filename+dims vector |
| `rsf_read.c` | MEX: read samples from an RSF file into a pre-allocated MATLAB double buffer |
| `rsf_read_all.c` | MEX: read an entire RSF file in one call, returning data and optional header metadata |
| `rsf_read_header.c` | MEX: read RSF header only (sizes, deltas, origins, labels, units) without touching the binary |
| `rsf_write.c` | MEX: write a MATLAB double array to an RSF file, with optional `'same'` append mode |
| `rsf_write_all.c` | MEX: write a complete RSF file (header + data) in one call with full metadata support |
| `rsf_par.c` | MEX: read a scalar parameter (int, float, or bool) from an RSF header by name |
| `rsf_dim.c` | MEX: return the axis-length vector `[n1; n2; ...]` for an RSF file |
| `m8r.c` | MEX: run any Madagascar `sf*` binary as a filter on MATLAB data via temp RSF files |
| `test/clip.m` | Example MATLAB sf-program: clips data to `[-clip, clip]` trace by trace using `rsf_dim`, `rsf_create`, `rsf_read`/`rsf_write` with `'same'` |

---

## Shared conventions

All Madagascar sf-programs — regardless of language — follow the same conventions for
file naming, self-documentation comments, parameter style, error handling, and testing.
See [`skills/authoring-sf-programs/`](../authoring-sf-programs/SKILL.md) for the
full shared reference.
