---
name: authoring-sf-programs-f90
description: Use when authoring a Madagascar sf* main program in Fortran 90.
---

## When to use

Load this skill when writing a new `sf<name>` main program in **Fortran 90**. F90 is preferred when:

- The algorithm already exists as Fortran legacy or numerical code and a rewrite would be risky.
- The researcher team works primarily in Fortran and a C translation would lose maintainability.
- Dense array operations benefit from Fortran's natural multi-dimensional array semantics.

The source file naming convention is `M<name>.f90`, placed inside `user/<youruser>/`. The `M` prefix tells the SCons build system that this is a main program; the installed binary will be named `sf<name>`.

**Dependency note:** The F90 API requires a Fortran 90 compiler. `gfortran` is the reference compiler; `pgf90` and `ifort` are also supported (see `SConstruct` logic). The Madagascar `configure` step must detect an F90 compiler and set `API=f90`; if detection fails the `librsff90` library is not built and F90 programs cannot be compiled.

For language-agnostic conventions (self-documentation comments, parameter style, build integration, error handling, testing), read the shared skill first:

- [skills/authoring-sf-programs/](../authoring-sf-programs/SKILL.md)

---

## Skeleton

Minimal F90 program structure, derived from `api/f90/Testfile.f90`:

```fortran
! sf<name> — one-line description
!
! Usage: sf<name> < inp.rsf > out.rsf [param=value ...]
program Mname
  use rsf

  implicit none
  type(file)                       :: in, out
  integer                          :: n1, n2, i2
  real                             :: d1, o1
  real, dimension(:), allocatable  :: trace

  ! 1. Initialise RSF and parse command-line parameters
  call sf_init()

  ! 2. Open standard input / output files
  in  = rsf_input()          ! default tag "in"
  out = rsf_output("out")    ! default tag "out"

  ! 3. Read axis metadata from input history
  call from_par(in, "n1", n1)          ! fast axis length (required)
  call from_par(in, "d1", d1, 1.0)    ! sampling (optional, default 1)
  call from_par(in, "o1", o1, 0.0)    ! origin  (optional, default 0)
  call from_par(in, "n2", n2, 1)      ! slow axis (optional, default 1)

  ! 4. Read user parameter
  ! call from_par("key", val)  -- reads from command line; aborts if missing
  ! call from_par("key", val, default)  -- uses default if absent

  ! 5. Propagate / update output history
  call to_par(out, "n1", n1)
  call to_par(out, "n2", n2)

  ! 6. Allocate working array
  allocate(trace(n1))

  ! 7. Trace loop: read → process → write
  do i2 = 1, n2
     call rsf_read(in,  trace)
     ! ... process trace ...
     call rsf_write(out, trace)
  end do

end program Mname
```

Key structural points:

- `use rsf` is the only module import needed; it re-exports the entire public API.
- `implicit none` is mandatory — the module itself uses it, and omitting it in the main program will cause silent type bugs.
- `rsf_input()` / `rsf_output()` accept an optional `tag` argument (default `"in"` and `"out"` respectively) that corresponds to the pipe or file descriptor label used on the command line.
- `from_par` and `to_par` are generic interfaces; the compiler dispatches on the type of the `value` argument (`integer`, `real`, `character`, `logical`, arrays).
- Allocate arrays *after* reading axis dimensions, never before.

---

## API cheat sheet

All entries are sourced from `api/f90/rsf.f90` and `api/f90/fortran.c`. "Subroutine" means called with `call`; "function" means used in an expression or assignment.

| Purpose | F90 call | Notes |
|---------|----------|-------|
| **Initialise** | `call sf_init()` | Must be first; parses argv, sets up pipes |
| **Open input** | `in = rsf_input([tag])` | Function returns `type(file)`; default tag `"in"` |
| **Open output** | `out = rsf_output([tag])` | Function returns `type(file)`; default tag `"out"` |
| **Get hist int** | `call from_par(f, "n1", ival [,default])` | Reads integer from file header; aborts if missing and no default |
| **Get hist real** | `call from_par(f, "d1", rval [,default])` | Reads real from file header |
| **Get hist string** | `call from_par(f, "key", sval [,default])` | String result; max `FSTRLEN=256` chars |
| **Get hist int array** | `call from_par(f, "key", ivec [,default])` | `ivec` is `integer, dimension(:)` |
| **Get param int** | `call from_par("key", ival [,default])` | No `file` arg — reads from command line |
| **Get param real** | `call from_par("key", rval [,default])` | No `file` arg — reads from command line |
| **Get param bool** | `call from_par("key", lval [,default])` | `lval` is `logical` |
| **Get param real array** | `call from_par("key", rvec [,default])` | `rvec` is `real, dimension(:)` |
| **Get param string** | `call from_par("key", sval [,default])` | String from command line |
| **Put hist int** | `call to_par(out, "n1", ival)` | Writes integer key into output header |
| **Put hist real** | `call to_par(out, "d1", rval)` | Writes real key into output header |
| **Put hist string** | `call to_par(out, "key", sval)` | Writes string key into output header |
| **Put hist int array** | `call to_par(out, "key", ivec)` | Writes integer array |
| **Read float data** | `call rsf_read(in, array)` | Dispatches on rank (1-D to 5-D); also accepts explicit `n` |
| **Write float data** | `call rsf_write(out, array)` | Same dispatch rules as `rsf_read` |
| **Read complex data** | `call rsf_read(in, carray)` | `carray` is `complex, dimension(:...)` |
| **Write complex data** | `call rsf_write(out, carray)` | `carray` is `complex, dimension(:...)` |
| **Get/put full axis** | `call iaxa(f, ax, i)` / `call oaxa(f, ax, i)` | Reads/writes `type(axa)` struct (n,o,d,label,unit) for axis `i` |
| **Print axis** | `call raxa(ax)` | Writes axis summary to `stderr` |
| **Get file type** | `t = gettype(f)` | Returns `sf_float=3`, `sf_int=2`, `sf_complex=4`, etc. |
| **Set file type** | `call settype(f, t)` | Pass type constant, e.g. `sf_int` |
| **File size** | `s = filesize(f [,dim])` | Without `dim`: total elements; with `dim`: elements from axis `dim` onward |
| **File dimensions** | `nd = dimension(f, n)` | Fills integer array `n`; returns number of dims |
| **Seek** | `call rsf_seek(f, offset, whence)` | `whence`: `sf_seek_set=0`, `sf_seek_cur=1`, `sf_seek_end=2` |
| **Flush** | `call rsf_fileflush(out [,src])` | Flush output; optionally copy header from `src` |
| **Fatal error** | `call sf_error("message")` | Prints to stderr and aborts |

Type constants declared in `rsf.f90`: `sf_uchar=0`, `sf_char=1`, `sf_int=2`, `sf_float=3`, `sf_complex=4`, `sf_short=5`.

Seek constants: `sf_seek_set=0`, `sf_seek_cur=1`, `sf_seek_end=2`.

---

## Build integration

The F90 API library is built by `api/f90/SConstruct`:

1. It checks `'f90' in env.get('API',[]) and 'F90' in env` — both must be true, which requires that `configure` detected a Fortran 90 compiler.
2. It compiles `fortran.c` (the C-to-Fortran glue) and `rsf.f90` into `librsff90.a`.
3. It installs `librsff90.a` → `lib/` and the compiled module file `rsf.mod` → `include/`.
4. Compiler-specific flags are set automatically:
   - `gfortran` / `gfc`: adds `-DGFORTRAN` (selects `_gfortran_iargc` / `_gfortran_getarg_i4` symbols in `fortran.c`).
   - `pgf90`: extra include path only.
   - `ifort`: adds `-DINTEL_COMPILER`.

For user F90 programs (`user/<youruser>/M<name>.f90`), follow the same pattern as for C/C++ programs but use the `f90` target helpers described in the shared skill. The SCons infrastructure links against `librsff90` automatically when it detects a `.f90` source file with an `M` prefix. See the shared skill ([skills/authoring-sf-programs/](../authoring-sf-programs/SKILL.md)) for the `SConstruct` snippet and `UserSconsTargets.f90` details.

---

## Pointers

Files in `api/f90/` (as of this writing):

| File | Description |
|------|-------------|
| `rsf.f90` | The RSF Fortran-90 module; defines `module RSF`, `type(file)`, `type(axa)`, and all public generic interfaces (`from_par`, `to_par`, `rsf_read`, `rsf_write`, etc.) |
| `fortran.c` | C glue layer; wraps C API functions (`sf_init`, `sf_histint`, `sf_floatread`, etc.) into Fortran-callable entry points using `cfortran.h` macros |
| `cfortran.h` | Portable C-to-Fortran calling-convention header (third-party); handles name mangling for gfortran, pgf90, ifort, and others via `FCALLSC*` / `CCALLSF*` macros |
| `SConstruct` | SCons build script; detects F90 compiler, compiles `librsff90.a`, installs library and module, and builds the two test programs |
| `Testfile.f90` | Template/integration test for file I/O (`rsf_input`, `rsf_output`, `from_par`, `rsf_read`, `rsf_write`) — the canonical minimal program skeleton |
| `Testgetpar.f90` | Template/integration test for parameter parsing (`from_par` on command-line args, bool, real arrays) and `sf_error` |
| `test` | Directory with additional regression-test helpers |

---

## Shared conventions

All conventions below are defined in the shared skill and apply unchanged to F90 programs:

- Self-documentation comment block at the top of the file (processed by `doc/`)
- Parameter naming and ordering rules
- Error message style (`sf_error` with lowercase message, no trailing period)
- Testing via `sftour` / `SConstruct` `Test` targets
- Build file location (`user/<youruser>/SConstruct` or top-level `SConstruct`)

See [skills/authoring-sf-programs/](../authoring-sf-programs/SKILL.md) for the full reference.
