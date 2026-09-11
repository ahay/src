---
name: authoring-sf-programs-cpp
description: Use when authoring a Madagascar sf* main program in C++.
---

## When to use

Load this skill when you are writing a new `sf<name>` main program in C++. Choose C++ over C when:

- The algorithm is naturally expressed using **templates** (generic numerical kernels, policy-based designs, expression-template math libraries).
- You are integrating with an existing **C++ library ecosystem** (Eigen, Boost, RVL/RVLCOM operator framework, or similar).
- You want **object-oriented encapsulation** of state across multiple passes or operator classes.
- You are porting **existing scientific C++ code** and the translation cost to C is not justified.

The source file for a C++ `sf<name>` program is named `M<name>.cc` and lives in `user/<youruser>/`. The installed binary is `sf<name>` (the `M` is dropped, `sf` is prepended) — identical convention to C.

This skill is C++-specific. For language-agnostic conventions — file naming, the self-documentation comment format, parameter conventions, error handling, and testing — see the companion skill:

- [`../authoring-sf-programs/SKILL.md`](../authoring-sf-programs/SKILL.md) — shared conventions (load this too)

---

## Skeleton

The structure below follows `api/c++/Testfile.cc` exactly. Copy it verbatim and extend from here.

```cpp
// One-sentence description of what this program does.

#include <valarray>
#include <rsf.hh>

int main(int argc, char* argv[])
{
    sf_init(argc, argv);

    iRSF par(0);   // parameter object (0 = command-line only, not a file)
    iRSF in;       // default: opens "in" (stdin)
    oRSF out;      // default: opens "out" (stdout)

    // Read axis metadata from the input header
    int n1;
    float d1, o1;
    in.get("n1", n1);
    in.get("d1", d1);
    in.get("o1", o1);

    // Read a command-line parameter; supply default with three-arg form
    int n2;
    par.get("n2", n2, 1);
    /* n2 — number of output slices */

    // Check data type if it matters
    if (in.type() != SF_FLOAT)
        sf_error("Need float input.");

    // Allocate a trace buffer using std::valarray
    std::valarray<float> trace(n1);

    // Write axis metadata to the output header before the data loop
    out.put("n1", n1);
    out.put("d1", d1);
    out.put("o1", o1);
    out.put("n2", n2);

    // Main I/O loop: read one trace, process, write
    for (int i2 = 0; i2 < n2; i2++) {
        in >> trace;
        // ... process trace ...
        out << trace;
    }

    exit(0);
}
```

Key points:

- `sf_init(argc, argv)` must be the first call.
- `iRSF par(0)` opens the parameter object in command-line-only mode (not a file). Use `iRSF par` (no argument) only when `par` is also a data file.
- `in.get("n1", n1)` reads axis metadata (header key `n1`). No default — aborts if absent.
- `par.get("key", var, default)` reads a command-line parameter with a default.
- Allocate buffers with `std::valarray<T>` for automatic memory management.
- Write all output header keys with `out.put(...)` **before** the data loop.
- Use `>>` to read and `<<` to write `std::valarray` buffers.
- `sf_error("message")` (the C function, accessible via `rsf.hh`) prints to stderr and exits.

---

## API cheat sheet

All calls below are derived directly from `api/c++/rsf.hh`.

| Operation | C++ call |
|-----------|----------|
| Initialize | `sf_init(argc, argv);` |
| Open default input | `iRSF in;` — opens `"in"` (stdin) |
| Open named input | `iRSF vel("vel");` — opens file passed as `vel=` |
| Open parameter object | `iRSF par(0);` — command-line params, no file |
| Open default output | `oRSF out;` — opens `"out"` (stdout) |
| Open named output | `oRSF wt("weight");` — opens file passed as `weight=` |
| Read axis int (required) | `in.get("n1", n1);` |
| Read axis float (required) | `in.get("d1", d1);` |
| Read axis string | `std::string label; in.get("label1", label);` |
| Write axis int | `out.put("n1", n1);` |
| Write axis float | `out.put("d1", d1);` |
| Write axis string | `out.put("label1", "Time");` |
| Read int param (required) | `par.get("niter", niter);` |
| Read int param (with default) | `par.get("niter", niter, 100);` |
| Read float param | `par.get("eps", eps, 0.01f);` |
| Read bool param | `par.get("adj", adj, false);` |
| Read string param | `std::string mode; par.get("mode", mode, std::string("exact"));` |
| Read data (valarray) | `in >> trace;` where `trace` is `std::valarray<float>` |
| Write data (valarray) | `out << trace;` |
| Read scalar | `float v; in >> v;` |
| Write scalar | `float v = 1.f; out << v;` |
| Set output data type | `out.type(SF_INT);` |
| Check input data type | `if (in.type() != SF_FLOAT) sf_error("need float");` |
| Query total file size | `int total = in.size(0);` — product of all axes |
| Query size along axis k | `int nk = in.size(k);` — size of axis k (1-based) |
| Error handler | `sf_error("msg: %d", val);` — stderr + exit |

Notes on `put` overloads in `oRSF`: there are three overloads — `put(name, int)`, `put(name, float)`, and `put(name, const char*)`. There is **no** `put(name, float, size, array)` overload for float arrays (that overload is commented out in `rsf.hh`); use the int-array form `put(name, size, int_array)` only.

---

## Build integration

### What `api/c++/SConstruct` does

`api/c++/SConstruct` compiles `rsf.cc` and `cub.cc` into a static library named `rsf++` (file: `librsf++.a`, installed to `lib/`). The relevant line:

```python
lib = env.StaticLibrary('rsf++', ccfiles, CCFLAGS='')
env.Install('../../lib', lib)
env.Install('../../include', hhfiles)   # installs rsf.hh and cub.hh
```

It also prepends `../../include` to `CPPPATH` and `../../lib` to `LIBPATH`, and links against `librsf` (the C core). Test programs (`Testfile.x`, `Testgetpar.x`) are built in-place for local verification.

### For your user-directory C++ program

Use `HuiSconsTargets` in your `user/<youruser>/SConstruct`, **not** `UserSconsTargets` (the latter covers only `.c`, `.py`, `.f90`, `.jl`). The `HuiSconsTargets` helper exposes a `.cc` attribute:

```python
import sys, os
sys.path.append('../../framework')
import bldutil

targets = bldutil.HuiSconsTargets()
targets.cc = 'myprogram anotherprogram'   # base names without M prefix or .cc
targets.build_all(env, glob_build, srcroot, bindir, libdir, pkgdir)
```

This compiles `Mmyprogram.cc` → `sfmyprogram` and links it against both `librsf++` and `librsf`.

If your user directory does not yet have a `SConstruct`, copy the one from a nearby C++ user directory (e.g., `user/pyang/SConstruct` or `user/chenyk/SConstruct`) and adjust the program list.

### Linking flags (manual reference)

If you ever need to link manually outside SCons:

```bash
g++ -I$RSFROOT/include Mmyprogram.cc -L$RSFROOT/lib -lrsf++ -lrsf -o sfmyprogram
```

The C++ library is `-lrsf++` (from `librsf++.a`); the C core is `-lrsf`. Both must appear.

---

## Pointers to existing templates

All files in `api/c++/`:

- `Testfile.cc` — minimal I/O: read an `SF_INT` trace from stdin, write it five times to stdout. The simplest possible complete program. Start here.
- `Testgetpar.cc` — parameter parsing: demonstrates `par.get` for int, float, bool, and array variants with and without defaults.
- `rsf.hh` — the C++ API public header: `iRSF` and `oRSF` class declarations, all `get`/`put`/`>>` / `<<` overloads.
- `rsf.cc` — implementation of `iRSF` and `oRSF`; wraps `sf_input`, `sf_output`, `sf_histint`, `sf_histfloat`, `sf_getint`, `sf_getfloat`, `sf_getbool`, `sf_getstring`, `sf_floatread`, `sf_floatwrite`, etc.
- `cub.hh` — higher-level `CUB` class: manages an `sf_axis*` array for multi-dimensional cubes; exposes `headin()`, `headou()`, `clone()`, `getax(int)`, `putax(int, sf_axis)`, `setup(int kd)`, and typed `>>` / `<<` operators for float, int, short, char, sf_complex, and `std::complex<float>`.
- `cub.cc` — implementation of `CUB`; prefer `CUB` over raw `iRSF`/`oRSF` when you need per-axis `sf_axis` structs (origin, delta, label, unit).
- `SConstruct` — SCons build script: compiles `rsf.cc`+`cub.cc` into `librsf++.a`, installs `rsf.hh`+`cub.hh` to `include/`, builds test programs.
- `test/` — directory of additional test/regression scripts for the C++ API.

---

## Shared conventions

The self-documentation comment for a C++ program is a **single `//` line** immediately before the first `#include`. This is what `framework/rsf/doc.py` scrapes (`comment['c++']` regex). Parameter descriptions are `// trailing comments` on the same line as each `par.get(...)` call. File naming (`M<name>.cc`), parameter style (`key=value`), error handling (`sf_error`), and test patterns all follow the rules set out in the shared skill. For full details on all of these — including how `sfdoc` output is generated, how to write regression flows, and how to handle optional library dependencies — see [`skills/authoring-sf-programs/`](../authoring-sf-programs/SKILL.md).
