---
name: authoring-sf-programs-chapel
description: Use when authoring a Madagascar sf* main program in Chapel.
---

## When to use

Load this skill when writing a new `sf<name>` main program in Chapel. Chapel is a
parallel, HPC-oriented language developed by Cray/HPE. It is niche within the
Madagascar ecosystem: the `api/chapel/m8r.chpl` module exists and two test programs
(`api/chapel/test/clip.chpl`, `api/chapel/test/afdm.chpl`) exist, but there are
likely **no** `M*.chpl` user programs anywhere in the tree. Prefer C or Python
unless you specifically need Chapel's `forall` data-parallelism or locale-based
distributed computing.

Main programs use the extension `.chpl` and the file naming convention `M<name>.chpl`
(e.g., `user/yourname/Msmooth.chpl` installs as `sfsmooth`).

This skill is Chapel-specific. For language-agnostic conventions (file naming,
self-doc format, parameter style, build integration) also load the companion:

- [`../authoring-sf-programs/SKILL.md`](../authoring-sf-programs/SKILL.md) — shared conventions

---

## Skeleton

Minimal correct Chapel program. The entry point is `proc main(args: [] string)` —
Chapel passes command-line arguments as a string array, which is forwarded directly
to `sf_init`.

```chapel
// One-sentence description of what this program does.
use m8r;

proc main(args: [] string)
{
    // Initialize Madagascar (must be first RSF call)
    sf_init(args);

    // Open I/O files
    var fin:  sf_file = sf_input("in");
    var fout: sf_file = sf_output("out");

    // Read header metadata from input
    var n1: int(32);
    if !sf_histint(fin, "n1", n1) then
        sf_error("No n1= in input");

    var n2 = sf_leftsize(fin, 1);  // number of traces

    // Read a command-line parameter (required)
    var clip: real(32);
    if !sf_getfloat("clip", clip) then
        sf_error("Need clip=");

    // Allocate trace buffer
    var trace: [0..n1-1] real(32);

    // Main loop
    for i2 in 0..n2-1 {
        sf_floatread(trace, n1, fin);
        // ... process trace ...
        sf_floatwrite(trace, n1, fout);
    }

    // Close files
    sf_fileclose(fin);
    sf_fileclose(fout);

    // Finalize Madagascar
    sf_close();
}
```

Key invariants:

- `sf_init(args)` must be the **first** RSF call; pass `args` as-is.
- Open all `sf_input` files before `sf_output` — the first output inherits
  dimensions from the first input.
- Call `sf_close()` at the end to flush output headers.
- Use `sf_error("msg")` (not Chapel's builtin `halt`) for fatal errors —
  `sf_error` writes to stderr and exits with a non-zero status, which is the
  Madagascar pipeline convention.
- Chapel arrays passed to read/write routines must match the declared element
  type precisely: `real(32)` for float, `int(32)` for int, etc.

For parallel inner loops, replace the body loop with a `forall`:

```chapel
forall (ix, iz) in {0..<nx, 0..<nz} {
    ud[ix, iz] = ...;
}
```

---

## API cheat sheet

All procedures below are exposed by `api/chapel/m8r.chpl` via `use m8r;`.
The `//w` annotation in the source means a Chapel wrapper exists (accepts Chapel
strings / arrays directly); `//nw` means the raw C extern is used directly.

### Init / close

| Procedure | Signature | Notes |
|-----------|-----------|-------|
| `sf_init` | `(args: [] string)` | Must be first RSF call |
| `sf_close` | `()` | Flush headers; call at program end |

### File open / close

| Procedure | Returns | Notes |
|-----------|---------|-------|
| `sf_input(tag: string)` | `sf_file` | `"in"` = stdin |
| `sf_output(tag: string)` | `sf_file` | `"out"` = stdout |
| `sf_fileclose(file: sf_file)` | `void` | Close a file handle |

### File metadata

| Procedure | Notes |
|-----------|-------|
| `sf_gettype(file)` | Returns `sf_datatype` constant (`SF_FLOAT`, `SF_INT`, `SF_COMPLEX`, …) |
| `sf_settype(file, type_arg)` | Override data type on output |
| `sf_getform(file)` / `sf_setform(file, form)` | `SF_NATIVE`, `SF_ASCII`, `SF_XDR` |
| `sf_filesize(file): c_int` | Total element count |
| `sf_leftsize(file, dim: int(32)): int(32)` | Elements from dim onwards (use for trace count) |
| `sf_filedims(file, n: [] int(32)): int` | Fill array with all axis sizes |

### History (header) get — read from input file

| Procedure | Type read |
|-----------|-----------|
| `sf_histint(file, key, ref par: int(32)): bool` | integer |
| `sf_histfloat(file, key, ref par: real(32)): bool` | float |
| `sf_histdouble(file, key, ref par: real(64)): bool` | double |
| `sf_histbool(file, key, ref par: bool): bool` | boolean |
| `sf_histstring(file, key): c_string` | string (returns c_string) |
| `sf_histlargeint(file, key, ref par: int(64)): bool` | 64-bit int |
| `sf_histints(file, key, par: [] int(32), n): bool` | int array |
| `sf_histfloats(file, key, par: [] real(32), n): bool` | float array |
| `sf_histbools(file, key, par: [] bool, n): bool` | bool array |

### Header put — write to output file

| Procedure | Notes |
|-----------|-------|
| `sf_putint(file, key, par: c_int)` | write integer header |
| `sf_putfloat(file, key, par: c_float)` | write float header |
| `sf_putlargeint(file, key, par: c_long)` | write 64-bit int header |
| `sf_putstring(file, key, par: c_string)` | write string header |
| `sf_putints(file, key, par: [] c_int, n: c_int)` | write int array header |
| `sf_putline(file, line: c_string)` | write raw header line |

### Command-line parameter get

| Procedure | Type |
|-----------|------|
| `sf_getint(key, ref par: int(32)): bool` | integer |
| `sf_getfloat(key, ref par: real(32)): bool` | float |
| `sf_getdouble(key, ref par: real(64)): bool` | double |
| `sf_getbool(key, ref par: bool): bool` | boolean |
| `sf_getlargeint(key, ref par: int(64)): bool` | 64-bit int |
| `sf_getstring(key): c_string` | string |
| `sf_getints(key, par: [] int(32), n): bool` | int array |
| `sf_getfloats(key, par: [] real(32), n): bool` | float array |
| `sf_getbools(key, par: [] bool, n): bool` | bool array |
| `sf_getstrings(key, par: [] string, n): bool` | string array (colon-separated on CLI) |

### Data read / write

| Procedure | Chapel array type |
|-----------|-------------------|
| `sf_floatread(arr, size, file)` | `[] real(32)` |
| `sf_floatwrite(arr, size, file)` | `[] real(32)` |
| `sf_intread(arr, size, file)` | `[] int(32)` |
| `sf_intwrite(arr, size, file)` | `[] int(32)` |
| `sf_complexread(arr, size, file)` | `[] complex(64)` |
| `sf_complexwrite(arr, size, file)` | `[] complex(64)` |
| `sf_shortread(arr, size, file)` | `[] int(16)` |
| `sf_shortwrite(arr, size, file)` | `[] int(16)` |
| `sf_charread(arr, size, file)` | `[] int(8)` |
| `sf_charwrite(arr, size, file)` | `[] int(8)` |
| `sf_uncharread(arr, size, file)` | `[] uint(8)` |
| `sf_uncharwrite(arr, size, file)` | `[] uint(8)` |

### Diagnostics

| Procedure | Notes |
|-----------|-------|
| `sf_error(args...?n)` | Variadic; concatenates args, writes to stderr, exits non-zero. Use instead of `halt`. |
| `sf_warning(args...?n)` | Variadic; prints warning to stderr, continues execution. |

---

## Build integration

The `api/chapel/SConstruct` installs `m8r.chpl` into `$LIBDIR` so that user
programs can find the module with `-M$LIBDIR`. It does **not** define a
`UserSconsTargets.chpl` attribute — there is no automatic Chapel analog to the
C / Python / Fortran 90 discovery in `bldutil.UserSconsTargets()`.

To build a Chapel program manually, define a custom `Builder` in your user
directory's `SConstruct`, modelled on `api/chapel/test/SConstruct`:

```python
import rsf.proj
proj = rsf.proj.Project()

chprsf = Builder(
    action='$CHPL $CHPLFLAGS -I$INCDIR -L$LIBDIR -M$LIBDIR -l$LIBS $SOURCES $OPT -o $TARGET'
)
proj.Append(BUILDERS={'chprsf': chprsf})

proj.chprsf(
    'Msmooth.exe',
    ['Msmooth.chpl'],
    CHPL=proj.get('CHPL_HOST_COMPILER'),
    INCDIR=proj.get('CPPPATH'),
    LIBDIR=proj.get('LIBPATH'),
    LIBS='rsf',
    OPT=''
)
proj.End()
```

The Chapel compiler is available as `proj.get('CHPL_HOST_COMPILER')` when the
build system detected `chpl` at configure time. If `CHPL_HOST_COMPILER` is
`None`, Chapel was not found and you need to install it separately.

The resulting binary is named `Msmooth.exe` (or whatever `TARGET` you supply);
you can rename/install it manually as `sfsmooth`. Automatic installation via
`scons install` is not wired up for Chapel user programs.

---

## Pointers

Every file in `api/chapel/`:

| File | Description |
|------|-------------|
| `api/chapel/m8r.chpl` | The Chapel module (`module m8r`) that wraps `rsf.h`; exposes all RSF C functions as Chapel `extern proc` declarations plus Chapel-idiomatic wrapper procs for init, I/O, error, and data read/write. |
| `api/chapel/SConstruct` | SCons build script; installs `m8r.chpl` into `$LIBDIR` so user programs can `use m8r` via `-M$LIBDIR`. |
| `api/chapel/test/clip.chpl` | Minimal working example: reads float RSF, clips values to a user-supplied threshold, writes output. Good starting template. |
| `api/chapel/test/afdm.chpl` | Full parallel example: 4th-order finite-difference acoustic wave modeling; uses `forall` for 2-D stencil parallelism. |
| `api/chapel/test/SConstruct` | Builds `clip.exe` and `afdm.exe` using the custom `chprsf` Builder; includes a `Flow` regression test for `clip`. |

---

## Shared conventions

All language-agnostic rules (file naming, self-documentation, parameter style,
error handling, testing, build integration) live in:

- [`../authoring-sf-programs/SKILL.md`](../authoring-sf-programs/SKILL.md)

Chapel-specific reminders that differ from C conventions:

- **Error handling**: use `sf_error("msg")` (from `m8r`), not Chapel's builtin
  `halt`. `sf_error` is the Madagascar-standard way to abort a pipeline stage;
  it writes to stderr and exits non-zero so that SCons `Flow()` detects failure.
- **Entry point**: `proc main(args: [] string)` — not `int main(argc, argv)`.
  Pass `args` directly to `sf_init(args)`.
- **Type precision**: Chapel distinguishes `real(32)` vs `real(64)` and
  `int(32)` vs `int(64)`. RSF floats are 32-bit; use `real(32)` arrays with
  `sf_floatread` / `sf_floatwrite`. Mismatching widths causes a type error at
  compile time.
- **Array indexing**: Chapel arrays default to 1-based; declare domains
  explicitly (`[0..n1-1]`) to match RSF's 0-based C convention when interoping
  with `c_ptrTo`.
- **Parallelism**: `forall` loops are the primary Chapel parallelism construct
  and map naturally to RSF's regular-grid data model. `sf_floatread` and
  `sf_floatwrite` are serial I/O; parallelize the compute between them.
- **No `UserSconsTargets.chpl`**: Chapel user programs need a hand-written
  `SConstruct`; see **Build integration** above.
