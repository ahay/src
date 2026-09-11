---
name: authoring-sf-programs-f77
description: Use when authoring a Madagascar sf* main program in Fortran 77.
---

## When to use

Use this skill for Madagascar processing programs in **Fortran 77 fixed-format**
(file extension `.f`, program named `M<Name>.f`). Typical scenarios: maintaining
legacy fixed-format code, porting existing F77 algorithms, or working with
research groups that standardize on `gfortran -std=legacy` or `ifort`.

**F77 vs F90 key differences:**

| Aspect | F77 | F90 |
|---|---|---|
| Source format | Fixed-format; code starts at column 7 | Free-format |
| Continuation | Any non-space/zero char in column 6 | `&` at end of line |
| File extension | `M<name>.f` | `M<name>.f90` |
| Comments | `C` or `*` in column 1 | `!` anywhere |
| Loops | `do 10 i=1,n` … `10 continue` | `do i=1,n` … `end do` |

For shared conventions (SConstruct rules, self-doc headers, naming) see
[skills/authoring-sf-programs/](../authoring-sf-programs/SKILL.md).

---

## Skeleton

Minimal F77 program derived from `api/f77/Testfile.f` and `api/f77/test/clip.f`.
Columns 1-6 are reserved; statement body starts at column 7; continuation
marker goes in column 6.

```fortran
C     Mclip.f  --  clip float data to +/- clip
      program Clipit
      implicit none
      integer*8 in, out
      integer*8 sf_input, sf_output, sf_leftsize
      logical sf_getfloat, sf_histint
      integer n1, n2, i1, i2
      real clip, trace(1000)

      call sf_init()
      in  = sf_input("in")
      out = sf_output("out")

      if (.not. sf_histint(in,"n1",n1))
     &    call sf_error("No n1= in input")
      if (.not. sf_getfloat("clip",clip))
     &    call sf_error("Need clip=")
      n2 = sf_leftsize(in,1)

      do 10 i2=1, n2
         call sf_floatread(trace, n1, in)
         do 20 i1=1, n1
            if (trace(i1) .gt.  clip) trace(i1) =  clip
            if (trace(i1) .lt. -clip) trace(i1) = -clip
 20      continue
         call sf_floatwrite(trace, n1, out)
 10   continue
      stop
      end
```

Column layout reminder:
```
123456789...
C     comment (C in col 1)
      statement (col 7+)
     &continuation (col 6)
 10   label (cols 1-5)
```

---

## API cheat sheet

All calls go through the `librsff` wrapper (`fortran.c` + `cfortran.h`).
Functions returning a value must be declared with the appropriate Fortran type.

**Init and file handles** — declare file handles as `integer*8`:

| Call | Returns | Purpose |
|---|---|---|
| `call sf_init()` | — | Initialize Madagascar (must be first) |
| `sf_input("in")` | `integer*8` | Open named input RSF file |
| `sf_output("out")` | `integer*8` | Open named output RSF file |
| `call sf_fileclose(f)` | — | Close file handle |
| `call sf_close()` | — | Finalize I/O |

**History reads** — declare as `logical`; `.false.` means key absent:

| Call | Purpose |
|---|---|
| `sf_histint(f,"n1",i)` | Read integer header key |
| `sf_histfloat(f,"d1",x)` | Read float header key |
| `sf_histbool(f,"key",b)` | Read logical header key |
| `sf_histstring(f,"label")` | Read string (function returning `character*`) |

**History writes** (`call` subroutines, no return value):

| Call | Purpose |
|---|---|
| `call sf_putint(f,"n2",i)` | Write integer to output header |
| `call sf_putfloat(f,"d2",x)` | Write float |
| `call sf_putstring(f,"label",str)` | Write string |

**Command-line param reads** — declare as `logical`; `.true.` if key found:

| Call | Purpose |
|---|---|
| `sf_getint("niter",i)` | Read integer param |
| `sf_getfloat("clip",x)` | Read float param |
| `sf_getfloats("d",xv,n)` | Read n floats |
| `sf_getbool("verb",b)` | Read yes/no param |
| `sf_getstring("tag")` | Read string (function) |

**Data I/O** (`call` subroutines):

| Call | Purpose |
|---|---|
| `call sf_floatread(buf, n, f)` | Read `n` floats from `f` |
| `call sf_floatwrite(buf, n, f)` | Write `n` floats to `f` |
| `call sf_intread(buf, n, f)` | Read `n` integers |
| `call sf_intwrite(buf, n, f)` | Write `n` integers |
| `call sf_complexread(buf, n, f)` | Read `n` complex values |
| `call sf_complexwrite(buf, n, f)` | Write `n` complex values |

**Utility:**

| Call | Returns | Purpose |
|---|---|---|
| `sf_leftsize(f, dim)` | `integer*8` | Elements past dimension `dim` |
| `sf_gettype(f)` | `integer` | Data type code (3 = float) |
| `call sf_error("msg")` | — | Print and abort |
| `call sf_warning("msg")` | — | Print, continue |

---

## Build integration

The F77 API builds only when `f77` appears in the SCons `API` list
(`api/f77/SConstruct`):

1. Detects compiler via `env.get('F77')`; sets `-DGFORTRAN` or
   `-DINTEL_COMPILER` accordingly.
2. Compiles `fortran.c` (with `cfortran.h`) into `fortran.o`.
3. Bundles into `librsff.a` (static library) and installs to `lib/`.
4. User programs link against both `librsff` and core `librsf`.

Typical user `SConstruct` snippet:

```python
env.Append(LIBS=['rsff', 'rsf'])
env.Program('sfclip', ['Mclip.f'])
# For strict fixed-format:
env.Replace(F77FLAGS='-std=legacy -ffixed-form')
```

`cfortran.h` handles name-mangling: `gfortran` (single underscore) vs
`f2c`-style (double underscore) via compile-time defines.

---

## Pointers

Files in `/Users/jgoai/m8r/src/api/f77/`:

| File | Description |
|---|---|
| `fortran.c` | C wrapper; `FCALLSC*` macros expose every `sf_*` function to F77 |
| `cfortran.h` | Third-party (Burkhard Burow/DESY) Fortran↔C name-mangling header |
| `Testfile.f` | Smoke-test: open in/out, read 100-float traces, write them back |
| `Testgetpar.f` | Exercises all param-read calls (`sf_getint`, `sf_getfloat`, etc.) |
| `SConstruct` | Builds `librsff.a`, installs to `lib/`, compiles test programs |
| `test/clip.f` | Complete real program: float clipping with header inspection |
| `test/SConstruct` | SCons build for the `test/` example programs |

---

## Shared conventions

General Madagascar conventions (self-documentation strings, parameter help,
SCons `Program` builder, installation targets) are in:

[skills/authoring-sf-programs/](../authoring-sf-programs/SKILL.md)

F77-specific reminders:

- Declare `sf_input`, `sf_output`, `sf_leftsize`, `sf_filesize` as `integer*8`
  — file handles and offsets are 64-bit on modern systems.
- Declare param-query functions (`sf_getint`, `sf_histfloat`, …) as `logical`
  to test them in `if (.not. …)` guards.
- Use `implicit none` — F77 implicit typing silently maps `i`-`n` to integer
  and the rest to real, hiding declaration bugs.
- Continuation: place the continuation character in **column 6** only.
