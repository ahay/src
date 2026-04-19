---
name: authoring-sf-programs-julia
description: Use when authoring a Madagascar sf* main program in Julia.
---

## When to use

Load this skill when writing a new `sf<name>` program in Julia. Julia is a good fit
when your collaborators already work in that ecosystem, when you want direct access to
the Julia numerical stack (FFTW.jl, LinearAlgebra, LoopVectorization, etc.) without a
foreign-function call, or when you need seamless FFI to C but prefer Julia syntax over
Python.

The source file is named `M<name>.jl` and lives in `user/<youruser>/`. The build
system handles the rest. Framework support for Julia programs exists via
`UserSconsTargets.jl` (see **Build integration**). There are no existing `M*.jl`
programs in the tree yet; use `api/julia/test/` as templates (especially `clip.jl`
for the core pattern and `afdm.jl` for a heavier numerical example).

The Julia API (`m8r.jl`) is a thin wrapper around `libdrsf` via `ccall`. All data
types map to their C equivalents: `Float32`, `Int32`, `ComplexF32`, `Int16`, `UInt8`.

See also: [Shared language-agnostic conventions](../authoring-sf-programs/SKILL.md)

---

## Skeleton

Modelled on `api/julia/test/clip.jl`:

```julia
#!/usr/bin/env julia
# M<name>.jl — one-sentence description of what this program does.
#
# Additional comments, references, and notes go here.
# (No structured doc-scraper for Julia yet; the comment is for human readers.)

import m8r

# --- Parameter retrieval (from command line or par= file) ---
clip  = m8r.getfloat("clip")          # required — exits if absent
verb  = m8r.getbool("verb",  false)   # optional — default false
niter = m8r.getint("niter",  100)     # optional — default 100

# --- Open standard input / output ---
inp = m8r.input("in")
out = m8r.output("out")

# --- Read header fields from inp ---
n1 = m8r.histint(inp, "n1")
n2 = m8r.leftsize(inp, 1)   # total samples / n1; works for any nd array

# --- Copy header to output, optionally modifying axes ---
m8r.putint(out,    "n1", n1)
m8r.putfloat(out,  "d1", m8r.histfloat(inp, "d1"))
m8r.putfloat(out,  "o1", m8r.histfloat(inp, "o1"))
m8r.putstring(out, "label1", m8r.histstring(inp, "label1"))
m8r.putstring(out, "unit1",  m8r.histstring(inp, "unit1"))

# --- Trace loop ---
trace = Array{Float32}(undef, n1)
for i2 in 1:n2
    m8r.floatread(trace, n1, inp)

    # ... process trace ...
    clamp!(trace, -clip, clip)

    m8r.floatwrite(trace, n1, out)
end
```

Key structural points:
- `import m8r` (not `using m8r`) keeps the namespace explicit; `using m8r` is also
  valid and exports `rsf_read` / `rsf_write`.
- `m8r.__init__()` is called automatically when the module loads; do not call
  `sf_init` manually.
- `m8r.leftsize(file, 1)` returns `total_samples / n1` — the product of all axes
  beyond the first. Use this to drive the outer loop without hard-coding dimensions.
- Auxiliary files: `m8r.input("vel")` opens the file named by the command-line
  argument `vel=<path>`.

---

## API cheat sheet

| Purpose | Julia call | Notes |
|---------|-----------|-------|
| **Init** | automatic on `import m8r` | `__init__` calls `sf_init` via `ccall` |
| **Open input** | `m8r.input("in")` | `"in"` = stdin; any other string = file path from that CLI param |
| **Open output** | `m8r.output("out")` | `"out"` = stdout; other strings = named output files |
| **Close file** | `m8r.close(file)` | wraps `sf_fileclose`; called automatically by `rsf_write` |
| **Hist int** | `m8r.histint(file, "n1")` | read `Int` from header; returns `0` if key absent |
| **Hist float** | `m8r.histfloat(file, "d1")` | read `Float32` from header |
| **Hist string** | `m8r.histstring(file, "label1")` | read string from header; `""` if absent |
| **Put int** | `m8r.putint(out, "n1", n)` | write integer header key |
| **Put float** | `m8r.putfloat(out, "d1", d)` | write float header key |
| **Put string** | `m8r.putstring(out, "label1", s)` | write string header key |
| **Get int param** | `m8r.getint("niter", 100)` | CLI/par param; default if absent |
| **Get float param** | `m8r.getfloat("clip")` | returns `0f0` if absent and no default given |
| **Get string param** | `m8r.getstring("vel", "")` | returns default string if absent |
| **Get bool param** | `m8r.getbool("verb", false)` | `"y"` → `true`, `"n"` → `false` |
| **Read floats** | `m8r.floatread(arr, n, inp)` | `arr` must be `Array{Float32,1}` of length `n` |
| **Write floats** | `m8r.floatwrite(arr, n, out)` | `arr` must be `Array{Float32,1}` |
| **Read ints** | `m8r.intread(arr, n, inp)` | `Array{Int32,1}` |
| **Write ints** | `m8r.intwrite(arr, n, out)` | `Array{Int32,1}` |
| **Read complex** | `m8r.complexread(arr, n, inp)` | `Array{ComplexF32,1}` |
| **Write complex** | `m8r.complexwrite(arr, n, out)` | `Array{ComplexF32,1}` |
| **Read short** | `m8r.shortread(arr, n, inp)` | `Array{Int16,1}` |
| **Write short** | `m8r.shortwrite(arr, n, out)` | `Array{Int16,1}` |
| **Left size** | `m8r.leftsize(file, dim)` | samples in dims `> dim`; `leftsize(f,0)` = total |
| **Set format** | `m8r.setformat(out, "complex")` | must call before first write for non-float output |
| **High-level read** | `dat, n, d, o, l, u = rsf_read(file_or_name)` | returns array + axis metadata |
| **High-level write** | `rsf_write(name, dat)` | writes array to named RSF file |
| **Pipe sf programs** | `sfspike(n1=10) \|> sfsmooth \|> rsf_read` | all installed sf* programs are exported |

---

## Build integration

`api/julia/SConstruct` installs `m8r.jl` into `$RSFROOT/lib/` at build time. The
framework-level dispatch for Julia programs happens in `UserSconsTargets.jl`
(located in the SCons framework tree). To add a Julia main program to your user
directory's build:

```python
# user/<youruser>/SConstruct
import sys, os
sys.path.append('../../framework')
import bldutil

# New-style (preferred)
targets = bldutil.UserSconsTargets()
targets.jl = 'myprogram'          # base name only — no M prefix, no .jl suffix
targets.build_all(env, glob_build, srcroot, bindir, libdir, pkgdir)
```

At configure time, Madagascar's build system checks for a Julia executable. If Julia
is not found, `targets.jl` programs are silently skipped; no error is raised. To
verify Julia was detected during configuration, inspect `$RSFROOT/include/config.py`
for a `JULIA` key, or check `env.get('JULIA')` inside the SConstruct.

Self-documentation for Julia programs is not yet scraped by `framework/rsf/doc.py`
(no `comment['jl']` regex entry); write a clear comment block at the top of your
`M<name>.jl` for human readers, and plan to add `sfdoc` support manually if needed.

---

## Pointers

Files in `api/julia/` and their purpose:

| File | Description |
|------|-------------|
| `api/julia/m8r.jl` | Full Julia API module — the only file you `import`; wraps all C-API entry points via `ccall` and also auto-exports every installed `sf*` binary as a Julia function |
| `api/julia/SConstruct` | SCons build file for the API; installs `m8r.jl` to `$RSFROOT/lib/` |
| `api/julia/test/clip.jl` | Minimal template: open in/out, read header, loop over traces with `floatread`/`floatwrite` — the closest thing to an `M*.jl` skeleton |
| `api/julia/test/afdm.jl` | Full numerical example: acoustic finite-difference modelling; shows multi-file input, header axis assembly, `putint`/`putfloat`/`putstring` output header setup, and a `@fastmath @inbounds` compute kernel |
| `api/julia/test/runtests.jl` | Comprehensive unit tests for every type (`uchar`, `char`, `int`, `float`, `complex`, `short`) covering `getint`/`getfloat`/`getstring`/`getbool`, `histint`/`histfloat`/`histstring`, `putint`/`putfloat`/`putstring`, read/write round-trips, and high-level `rsf_read`/`rsf_write` plus pipe-based `sf*` function calls |

---

## Julia-specific quirks

- **Type precision**: the API uses `Float32` and `Int32` throughout (matching the C
  API). Passing `Float64` arrays to `floatwrite` requires explicit conversion:
  `m8r.floatwrite(Float32.(vec(arr)), n, out)`. The high-level `rsf_write` does this
  conversion automatically.
- **Output type inheritance**: `m8r.output("out")` silently inherits the data type
  from the most recent `m8r.input("in")` call (C API behaviour). For non-float output
  (complex, int, short) call `m8r.setformat(out, "complex")` immediately after
  `m8r.output` and before any write. The `rsf_write(name, arr)` variant avoids this
  pitfall by writing a dummy pipe to force the correct type.
- **No `sf_error` in Julia**: there is no `m8r.error(...)` binding. Use `error("msg")`
  (Julia built-in) or `throw(...)` for fatal errors; the non-zero exit will propagate
  correctly through Madagascar pipelines.
- **Pipe composition**: all installed `sf*` binaries are available as Julia functions
  after `using m8r`. They accept an `RSFFile`, an array, or no argument, and return
  an `RSFFile` pointing to a temporary file. Chain them with `|>`.
- **RSFROOT required**: the module reads `ENV["RSFROOT"]` at load time. If the
  environment variable is not set, `m8r.RSFROOT` is `nothing` and `sf*` function
  generation is skipped; the low-level read/write API still works.

---

## Shared conventions

For file naming, self-documentation, parameter conventions, error handling, testing,
and SCons build patterns that apply to all languages, see:

[skills/authoring-sf-programs/](../authoring-sf-programs/SKILL.md)
