---
name: authoring-sf-programs-c
description: Use when authoring a new Madagascar sf* main program in C (the reference implementation — all other language APIs wrap this).
---

## When to use

Load this skill whenever you are writing a new `sf<name>` main program in C. C is the **reference implementation** for Madagascar: the Fortran 77, Fortran 90, C++, CUDA, and Python APIs all wrap or mirror the C API. Every claim in this skill is grounded in `api/c/rsf.h` (the amalgamated public header at `build/api/c/rsf.h`) and the test programs under `api/c/`.

This skill is C-specific. For language-agnostic conventions (file naming, self-doc format, parameter style, build integration) see the companion:

- `../authoring-sf-programs/SKILL.md` — shared conventions (load this too)

Typical placement: `user/<youruser>/M<name>.c`. The build system finds every `M*.c` in a user directory and compiles it into `sf<name>`.

---

## Skeleton

This is the minimal correct skeleton for a C main program. Copy it verbatim and extend from here.

```c
#include <rsf.h>

int main(int argc, char* argv[])
{
    sf_file in, out;
    int n1, n2;
    float d1, o1;
    float *trace;

    sf_init(argc, argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_histint  (in, "n1", &n1)) sf_error("No n1= in input");
    if (!sf_histint  (in, "n2", &n2)) n2 = 1;
    if (!sf_histfloat(in, "d1", &d1)) d1 = 1.0f;
    if (!sf_histfloat(in, "o1", &o1)) o1 = 0.0f;

    trace = sf_floatalloc(n1);

    for (int i2 = 0; i2 < n2; i2++) {
        sf_floatread(trace, n1, in);
        /* process trace here */
        sf_floatwrite(trace, n1, out);
    }

    free(trace);
    exit(0);
}
```

Key invariants:
- `sf_init` must be the first RSF call.
- Open all `sf_input` files before `sf_output` — the first output inherits dimensions from the first input by default.
- Always call `exit(0)` (not `return 0`) at the end; this is the Madagascar convention.

---

## Self-doc header

Every `M<name>.c` must begin with a comment block that the build system extracts to produce `sfdoc sf<name>` output. The format (sourced from `user/fomels/Mpick.c` in this tree):

```c
/* Automatic picking from semblance-like panels.

Takes: rect1=1 rect2=1 ...

rectN defines the size of the smoothing stencil in N-th dimension.

Theory in Appendix B of:
S. Fomel, 2009,
Velocity analysis using AB semblance: Geophysical Prospecting, v. 57, 311-321.
*/

/*
  Copyright (C) 2004 University of Texas at Austin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
```

Rules extracted from real programs in this tree:

1. **First comment** — brief one-line description on the opening line, optional extended prose below. This entire block becomes the program's documentation.
2. **Parameter inline comments** — each `sf_getint`/`sf_getfloat`/`sf_getbool` call is followed immediately by a `/* description */` comment on the next line. The doc extractor captures these and lists them as named parameters in `sfdoc`.
3. **Second comment** — the GPL copyright block. Keep it verbatim.
4. `#include <rsf.h>` comes after both comment blocks, before any other code.

Example parameter documentation from `user/fomels/Mpick.c`:

```c
if (!sf_getint("niter",&niter)) niter=100;
/* number of iterations */
if (!sf_getfloat("an",&an)) an=1.;
/* axes anisotropy */
if (!sf_getint("gate",&gate)) gate=3;
/* picking gate */
if (!sf_getbool("smooth",&smooth)) smooth=true;
/* if apply smoothing */
```

Each `/* ... */` directly following a `sf_get*` call is parsed as the parameter description; it appears under "parameters:" in `sfdoc sf<name>`.

---

## I/O API

Signatures from `build/api/c/rsf.h`.

### Initialization

```c
void sf_init(int argc, char *argv[]);
```
Must be the first RSF call. Parses the command line into an internal parameter table; opens stdin/stdout for piped RSF data.

### Opening files

```c
sf_file sf_input  (const char* tag);
sf_file sf_output (const char* tag);
```
`tag` is a logical name used on the command line: `sfmyprog in=data.rsf out=result.rsf`. Use `"in"` and `"out"` for the primary pipes (stdin / stdout). Additional named files (`"den"`, `"mask"`, etc.) are opened the same way. Output files inherit dimension headers from the first input by default — override with explicit `sf_put*` calls before the first write.

```c
sf_datatype sf_gettype(sf_file file);
```
Returns the data type (`SF_FLOAT`, `SF_COMPLEX`, `SF_INT`, etc.). Use to validate input type before reading.

```c
void sf_fileclose(sf_file file);
```
Closes an `sf_file` handle and flushes metadata. For secondary outputs (not stdin/stdout) it is good practice to call this explicitly; for primary stdin/stdout it is called implicitly on `exit`.

### Reading and writing float data

```c
void sf_floatread  (float* arr, size_t size, sf_file file);
void sf_floatwrite (float* arr, size_t size, sf_file file);
```
Read/write `size` consecutive floats. `arr` must be pre-allocated. These are the most common I/O calls.

### Reading and writing complex data

```c
void sf_complexread  (sf_complex* arr, size_t size, sf_file file);
void sf_complexwrite (sf_complex* arr, size_t size, sf_file file);
```
`sf_complex` is `float _Complex` (or the platform equivalent). Use when the input type is `SF_COMPLEX`.

---

## Header metadata API

RSF files carry key=value pairs in a text header. There are two directions:

- **Read from input header** — `sf_hist*` family. Returns `true` if the key exists.
- **Write to output header** — `sf_put*` family. Must be called before the first data write.

Signatures from `build/api/c/rsf.h`:

### Reading header values

```c
bool sf_histint    (sf_file file, const char* key, int* par);
bool sf_histfloat  (sf_file file, const char* key, float* par);
char* sf_histstring(sf_file file, const char* key);
```

`sf_histstring` returns a heap-allocated string (or `NULL`); test with `if (NULL != (...))`. Example:

```c
char *label;
if (NULL != (label = sf_histstring(scn, "label2")))
    sf_putstring(pik, "label", label);
```

### Writing header values

```c
void sf_putint    (sf_file file, const char* key, int par);
void sf_putfloat  (sf_file file, const char* key, float par);
void sf_putstring (sf_file file, const char* key, const char* par);
```

### Axis objects

Madagascar uses `sf_axis` structs to bundle `n`, `d`, `o`, `label`, and `unit` for each dimension. These simplify copying or transforming axes between files.

```c
sf_axis sf_iaxa(sf_file FF, int i);    /* read axis i (1-based) */
void    sf_oaxa(sf_file FF, const sf_axis AA, int i); /* write axis i */
```

Example — copy axis 1 from input to output unchanged, then set a new axis 2:

```c
sf_axis ax1, ax2;
ax1 = sf_iaxa(in, 1);
sf_oaxa(out, ax1, 1);
sf_putint(out, "n2", new_n2);
sf_putfloat(out, "d2", new_d2);
```

`sf_filedims` is a convenience wrapper that reads all `n` values at once:

```c
int sf_filedims(sf_file file, int *n);  /* returns number of dimensions */
```

---

## Parameter API

Command-line parameters are set by the caller as `key=value` pairs. Signatures from `build/api/c/rsf.h`:

### Scalar getters

```c
bool sf_getint      (const char* key, int* par);
bool sf_getfloat    (const char* key, float* par);
bool sf_getbool     (const char* key, bool* par);
char* sf_getstring  (const char* key);           /* returns NULL if absent */
bool sf_getlargeint (const char* key, off_t* par);
```

All scalar getters return `true` if the key was found (false otherwise). `sf_getstring` returns a pointer or `NULL`.

### Array getters

```c
bool sf_getints    (const char* key, int* par,   size_t n);
bool sf_getfloats  (const char* key, float* par, size_t n);
bool sf_getbools   (const char* key, bool* par,  size_t n);
```

Array getters read up to `n` comma-separated values (e.g. `rect=3,5,1`).

### Default-value pattern

The canonical Madagascar idiom — always provide a fallback default when the parameter is optional:

```c
if (!sf_getint("n", &n))     n = 100;
/* number of samples */
if (!sf_getfloat("eps", &eps)) eps = 0.01f;
/* regularization parameter */
if (!sf_getbool("verb", &verb)) verb = true;
/* verbosity flag */
```

The `/* ... */` on the line after each `sf_get*` call is parsed by the doc extractor and shown in `sfdoc sf<name>`.

For required parameters that have no sane default, call `sf_error` when absent:

```c
if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
```

---

## Memory API

Signatures from `build/api/c/rsf.h`. All allocators abort with a meaningful error message on allocation failure — never return `NULL`.

### 1-D allocators

```c
float*      sf_floatalloc  (size_t n);
sf_complex* sf_complexalloc(size_t n);
int*        sf_intalloc    (size_t n);
```

### 2-D and 3-D allocators

```c
float**      sf_floatalloc2 (size_t n1, size_t n2);
float***     sf_floatalloc3 (size_t n1, size_t n2, size_t n3);
sf_complex** sf_complexalloc2(size_t n1, size_t n2);
int**        sf_intalloc2   (size_t n1, size_t n2);
int***       sf_intalloc3   (size_t n1, size_t n2, size_t n3);
```

The multi-dimensional allocators return a pointer-to-pointer (standard C array-of-arrays), **but the underlying data is one contiguous block** — `out[0]` (for 2-D) or `out[0][0]` (for 3-D) points to it. This makes it safe to pass `arr[0]` to `sf_floatread`/`sf_floatwrite` as a flat buffer:

```c
float **scan = sf_floatalloc2(n1, n2);
sf_floatread(scan[0], n1 * n2, in);   /* reads entire 2-D block at once */
```

### Deallocation

Madagascar does **not** provide a custom deallocator. Use standard `free()` directly:

```c
free(trace);          /* 1-D */
free(scan[0]);        /* 2-D: free the contiguous data block */
free(scan);           /* then free the pointer array */
```

Most production Madagascar programs skip explicit `free` calls before `exit(0)` because the OS reclaims all memory on process exit. This is acceptable; add `free` calls when running under Valgrind or in regression test programs that check for leaks.

---

## Error API

Signatures from `build/api/c/rsf.h`:

```c
void sf_error  (const char *format, ...);
void sf_warning(const char *format, ...);
```

`sf_error` prints the formatted message to stderr, then calls `exit(1)`. Use it for unrecoverable conditions (bad input, missing required header key, wrong data type):

```c
if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
sf_error("Size mismatch [n%d]: %d != %d", j+1, m[j], n[j]);
```

`sf_warning` prints to stderr but **does not exit** — execution continues. Use it for progress messages and non-fatal conditions:

```c
sf_warning("cmp %d of %d;", i3+1, n3);  /* semicolon suppresses newline */
sf_warning(".");                          /* trailing dot flushes the line */
```

Both functions accept `printf`-style format strings and variadic arguments.

---

## Extra libraries

Most programs need only `#include <rsf.h>` and link against `librsf`. When a program needs helpers from `api/c/` (FFT, solvers, splines, etc.) or a user-directory helper module, the `SConstruct` in the user directory must declare the dependency.

### How `user/fomels/SConstruct` links extra libraries

The key lines from `user/fomels/SConstruct` (the real file in this tree):

```python
libs = [dynpre+'rsf'] + env.get('LIBS', [])

mains = Split(progs)
for prog in mains:
    sources = ['M' + prog]
    bldutil.depends(env, sources, 'M'+prog)
    env.StaticObject('M'+prog+'.c')
    prog = env.Program(prog, [x + '.o' for x in sources], LIBS=libs)
```

`bldutil.depends` scans the source file for `#include "foo.h"` lines and automatically adds `foo.o` to `sources`. This means that if `Mpick.c` includes `"dynprog.h"`, then `dynprog.c` is compiled and linked automatically — no manual addition to `sources` required.

For optional system libraries (FFTW, LAPACK, JPEG, TIFF), the pattern is:

```python
fftw = env.get('FFTW')
if fftw:
    env.Prepend(CPPDEFINES=['SF_HAS_FFTW'])

# LAPACK-dependent C++ programs
if lapack:
    libsxx = [dynpre+'rsf++', 'vecmatop']
    libsxx.extend(lapack)
    libsxx.extend(libs)
    prog = env.Program(prog, [x + '.cc' for x in sources], LIBS=libsxx)
```

For pure C programs that only need `librsf` helpers (the common case), the default `LIBS=libs` line is sufficient — add `#include "somehelper.h"` and `bldutil.depends` handles the rest.

---

## Worked references

All four files exist verbatim in this tree.

### `api/c/Testfile.c` — minimal I/O

Demonstrates the smallest possible main: `sf_init`, `sf_input`, `sf_output`, `sf_fileclose`, `exit`. No data is read or written. Useful to verify the build environment compiles cleanly.

### `api/c/Testgetpar.c` — parameter parsing

Shows `sf_getint`, `sf_getfloat`, `sf_getfloats`, `sf_getbool`, `sf_getbools`, and `sf_getstring` with assertion checks. The `par=Testgetpar.c` trick passes the source file itself as a par-file to supply default values — a technique used in regression tests.

### `api/c/Testfft.c` — FFT and `kiss_fft`

Uses `kiss_fft_next_fast_size` from `api/c/kiss_fft.h`, the portable FFT library bundled with Madagascar. Demonstrates including a sub-library header directly without going through `rsf.h`.

### `user/fomels/Mpick.c` — well-structured real-world program

Chosen because it demonstrates the full Madagascar C idiom in one file:

- Self-doc header with extended description and reference citation.
- `sf_init` / `sf_input` / `sf_output` pattern.
- `sf_histfloat` to read optional header metadata with defaults.
- `sf_getfloat`, `sf_getint`, `sf_getbool`, `sf_getstring` with inline doc comments.
- `sf_floatalloc2` (2-D allocation) and `sf_floatread`/`sf_floatwrite`.
- `sf_histstring` / `sf_putstring` to propagate axis labels.
- `sf_warning("cmp %d of %d;", ...)` progress reporting.
- `sf_unshiftdim` to reshape output dimensions.
- `exit(0)` at the end.

The only non-standard element is `#include "dynprog.h"` (a helper in the same directory), which `bldutil.depends` links automatically.

---

## Building and testing

### Build

Run `scons` inside `user/<youruser>/`:

```bash
cd /Users/jgoai/m8r/src/user/<youruser>
scons
```

The build system finds every `M*.c` in the directory automatically. No explicit registration in `SConstruct` is needed beyond listing the name in `progs`.

After `scons install` (or `make install` at the repo root), the binary is at:

```
/Users/jgoai/madagascar/bin/sf<name>
```

### Verify self-doc

```bash
sfdoc sf<name>
```

This confirms that the self-doc comment block was parsed correctly and that all `sf_get*` parameters appear with their inline descriptions.

### Verify basic I/O

```bash
echo "" | sfspike n1=10 | sf<name> > /dev/null
```

Or pipe through `sfin` to inspect the output header:

```bash
sfspike n1=100 n2=5 | sf<name> | sfin
```

### Regression test

For programs with non-trivial logic, write a `Testsf<name>.c` alongside `M<name>.c`. The pattern from `user/fomels/SConstruct`:

```python
for prog in Split('myalgorithm'):
    sources = ['Test' + prog, prog]
    bldutil.depends(env, sources, prog)
    sources = [x + '.o' for x in sources]
    env.Object('Test' + prog + '.c')
    env.Program(sources, PROGPREFIX='', PROGSUFFIX='.x', LIBS=libs)
```

The resulting `myalgorithm.x` binary can be run directly and is picked up by `scons` regression suites. See `api/c/Testfile.c` and `api/c/Testgetpar.c` for the testing style.

### Common pitfalls

- Forgetting `sf_init` → segfault on first `sf_input` call.
- Calling `sf_output` before the first `sf_input` → output header has no inherited dimensions.
- Reading `n1` from the command line instead of `sf_histint` → ignores the actual file dimensions.
- Using `malloc` instead of `sf_floatalloc` → no automatic out-of-memory error message, harder to debug.
- Omitting the inline `/* description */` after each `sf_get*` call → parameter missing from `sfdoc`.
