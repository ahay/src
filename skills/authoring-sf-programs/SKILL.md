---
name: authoring-sf-programs
description: Use when creating a new sf* main program — shared conventions that apply regardless of language (file naming, self-documentation, parameter style, build integration).
---

## When to use

Load this skill whenever you are writing a new `sf<name>` program, regardless of implementation language. It covers the language-agnostic layer: file naming, the self-documentation comment format, parameter conventions, error handling, testing, and build integration.

This skill is a **companion** to the language-specific skill `authoring-sf-programs-<lang>`. The per-language skill covers syntax — how to call `sf_init`, `sf_getbool`, `sf_error`, etc. in that language — while this skill establishes the shared contract that every sf-program must satisfy. Always load both:

1. `skills/authoring-sf-programs/SKILL.md` — you are here; shared conventions.
2. `skills/authoring-sf-programs-<lang>/SKILL.md` — language syntax and idioms.

---

## File naming

Every `sf<name>` main program lives in a source file named `M<name>.<ext>` inside a user directory `user/<youruser>/`. The `M` prefix is the signal to the build system that the file is a main program (not a library module). The installed binary is named `sf<name>` (the `M` is dropped, `sf` is prepended).

Extension to language mapping:

| Extension | Language |
|-----------|----------|
| `.c` | C |
| `.cc` | C++ |
| `.cu` | CUDA C |
| `.f` | Fortran 77 (framework support exists; no existing `M*.f` programs in this tree — see `api/f77/Test*.*` for templates) |
| `.f90` | Fortran 90 |
| `.py` | Python |
| `.jl` | Julia (framework support exists; no existing `M*.jl` programs in this tree — see `api/julia/Test*.*` for templates) |
| `.java` | Java (framework support exists; no existing `M*.java` programs in this tree — see `api/java/Test*.*` for templates) |
| `.m` | Matlab / Octave — disambiguated at configure time by which API is enabled (framework support exists; no existing `M*.m` programs in this tree — see `api/matlab/Test*.*` or `api/octave/Test*.*` for templates) |
| `.chpl` | Chapel |

Place your file at `user/<youruser>/M<name>.<ext>`. If the directory does not exist, create it along with a `SConstruct` (see **Build integration** below).

---

## Self-documentation

Every `sf<name>` program must begin with a structured comment block. The build system scrapes this block at compile/install time using `framework/rsf/doc.py` (specifically the `getprog()` function) and turns it into `sfdoc sf<name>` output.

### Comment format (C)

For C, the scraper uses this regex (from `framework/rsf/doc.py`, `comment['c']` regex):

```
comment['c'] = re.compile(r'\/\*(?P<comment>(?:[^*]+|\*[^/])+)\*\/')
```

That is: the **first** `/* ... */` block in the file. The format inside that block is:

- **Line 1**: one-sentence synopsis (becomes the `DESCRIPTION` section of `sfdoc`).
- **Remaining lines**: multi-line comments block (becomes `COMMENTS` in `sfdoc`). If the comments begin with `Takes: ...`, that suffix is appended to the synopsis line.

Real example from `user/fomels/Mpick.c` (which produces `sfdoc sfpick`):

```c
/* Automatic picking from semblance-like panels.

Takes: rect1=1 rect2=1 ...

rectN defines the size of the smoothing stencil in N-th dimension.

Theory in Appendix B of:
S. Fomel, 2009,
Velocity analysis using AB semblance: Geophysical Prospecting, v. 57, 311-321.
Reproducible version in RSFSRC/book/tccs/avo
http://ahay.org/RSF/book/tccs/avo/paper_html/

August 2012 program of the month:
http://ahay.org/blog/2012/08/01/program-of-the-month-sfpick/
*/
```

The first line `Automatic picking from semblance-like panels.` becomes the description. Everything else becomes the COMMENTS body.

For comparison, `system/main/transp.c` (the `sftransp` program):

```c
/* Transpose two axes in a dataset. 

If you get a "Cannot allocate memory" error, give the program a
memsize=1 command-line parameter to force out-of-core operation.
*/
```

The GPL license block that follows is a **separate** `/* ... */` comment and is not scraped.

### Parameter documentation (C)

Parameter documentation is harvested from `sf_get*` call sites. The scraper (`doc.py`, `param['c']` regex) matches:

```
if (!sf_get<type>("name",&var)) var=default;
/* description */
```

or, for required parameters (no default):

```
if (!sf_get<type>("name",&var)) sf_error("...");
```

The inline `/* ... */` comment immediately after the `sf_get` call becomes the parameter description in `sfdoc`. Example:

```c
if (!sf_getfloat("vel0",&vel0)) vel0=o2;
/* surface velocity */

if (!sf_getbool("smooth",&smooth)) smooth=true;
/* if apply smoothing */

if (!sf_getint("niter",&niter)) niter=100;
/* number of iterations */
```

For multi-valued (array) parameters, `sf_getints`, `sf_getfloats`, `sf_getbools`, `sf_getstrings` are matched by `params['c']` (note plural).

### Comment format in other languages

Each language has its own `comment[lang]`, `param[lang]`, and related regexes in `doc.py`:

- **Python**: first string literal (triple-quote or single-quote) in the file; parameters from `par.bool/int/float/string(...)` calls with an inline `# description` comment.
- **C++**: a single `// comment` line at the very top of the file (before any `#include`); parameters from `par.get(...)` calls with a trailing `// desc` comment.
- **Fortran 90**: `! comment` blocks; parameters from `from_par(...)` with `! desc`.
- **Chapel**: `// comment` lines; parameters from `config const/var name: type = default; // desc`.

The per-language skill has the precise syntax for each.

---

## Parameter conventions

These conventions are shared across all languages.

**Booleans** — always `y`/`n` on the command line. In C the type is `bool` and `sf_getbool` is used; the scraper normalizes `true`/`false` to `y`/`n` in the documentation. Users type `smooth=y` or `smooth=n`.

**Numeric defaults** — supply a reasonable default whenever the parameter is optional. When there is no sensible default, call `sf_error` (or its equivalent) if the parameter is absent. The sfdoc output shows the default in the synopsis.

**File-valued parameters** — auxiliary input/output files are opened with `sf_input("tag")` / `sf_output("tag")` where `tag` is anything other than `"in"` or `"out"`. The tag becomes the parameter name on the command line (`tag=filename.rsf`). The scraper detects these automatically from `inpout['c']` matches and adds them to the synopsis.

**Comma-separated lists** — use `sf_getints` / `sf_getfloats` / `sf_getbools` (plural) for array-valued parameters. The size argument is part of the scraper match and appears in the documentation.

**Required vs. defaulted** — if omission is fatal, write `if (!sf_getint("n",&n)) sf_error("n= required");`. If a default is acceptable, write `if (!sf_getint("n",&n)) n=100;`. The doc system distinguishes these by whether a default value appears in the scraped call.

Cross-reference: `skills/using-sf-programs/SKILL.md` explains how users supply parameters (`key=value` on the command line or via pipe from `sfdd` / a `.par` file).

---

## Error handling

**C**: `sf_error("format string: %d", value)` — prints the message to stderr and exits with a non-zero status. It is declared in `rsf.h` and defined in the core library.

**Python**: use `sys.stderr.write('msg\n'); sys.exit(1)` or `raise SystemExit('msg')`. (The RSF Python API does not provide a dedicated error helper — `sf_error` is a C API.)

**Fortran 90**: `call sf_error("msg")` (from `rsf_f90` module).

**C++**: `throw std::runtime_error("msg")` or `sf_error("msg")` (the C function is accessible from C++ via the C++ API headers).

**Chapel**: `sf_error("msg")` via the `RSF` module wrapper (see `api/chapel/m8r.chpl`).

**Pipeline behavior**: Madagascar programs communicate via Unix pipes (`stdin` → `stdout` in RSF format). When a program exits non-zero mid-pipe, the write end of the pipe closes. Downstream programs reading from stdin receive EOF and should terminate cleanly. If you run pipelines from a shell with `set -o pipefail`, any non-zero exit in the pipeline causes the entire pipeline to fail. In SCons `Flow()` commands Madagascar checks the exit status of each program automatically.

---

## Testing

**Low-level API tests** — `api/c/Test*.c` files test the C library itself (e.g., `api/c/Testfile.c`). These are compiled and linked by the `api/c/SConstruct` into `.x` test executables; they do not use the `M<name>.c` naming convention.

**Program-level tests** — in each user directory's `SConstruct`, test targets are SCons `Flow()` or `Command()` nodes. Run `scons` in `user/<youruser>/` to build and execute Flow/Command targets (including any test flows). The test typically pipes synthetic data through the new program and checks the output.

**Minimum test for a new program**: confirm that `sfdoc sf<name>` prints the correct description and parameter list. This validates that the self-doc comment was scraped correctly and that the program was installed.

**Regression test** (recommended): add a `Flow()` in your `SConstruct` that synthesizes known input with `sfspike` or `sfmath`, runs your program, and compares to a reference result with `sfattr` or a diff. This becomes part of the CI run.

CI runs the full test suite via `scons` from the build root; individual user directories can be tested with `scons` locally in `user/<youruser>/`.

---

## Build integration

**Standard case — no action required.** For C programs, the old-style `user/<author>/SConstruct` lists program base names in a `progs` string and loops:

```python
mains = Split(progs)
for prog in mains:
    sources = ['M' + prog]
    bldutil.depends(env, sources, 'M' + prog)
    env.StaticObject('M' + prog + '.c')
    prog = env.Program(prog, [x + '.o' for x in sources], LIBS=libs)
```

The new-style uses `bldutil.UserSconsTargets()`:

```python
targets = bldutil.UserSconsTargets()
targets.c = 'myprogram anotherprogram'
targets.py = 'myscript'
targets.build_all(env, glob_build, srcroot, bindir, libdir, pkgdir)
```

In both styles, you add the base name of your program (without `M` and without the extension) to the relevant list. The `SConstruct` handles compiling, linking against `librsf`, and installing to `bindir`.

For C++ (`.cc`) or CUDA (`.cu`) programs, use `HuiSconsTargets` instead — it exposes `.cc` and `.cu` attributes. `UserSconsTargets` covers only C (`.c`), Fortran-90 (`.f90`), Python (`.py`), and Julia (`.jl`).

**Self-doc** is generated by `env.Doc(prog, 'M' + prog)` (or `env.Doc(prog, 'M'+prog+'.py', lang='python')` for Python) and merged into the user's doc file. This depends on `framework/rsf/doc.py` and runs automatically as part of the build.

**Extra libraries** — if your program needs FFTW, for example, the `user/fomels/SConstruct` checks:

```python
fftw = env.get('FFTW')
if fftw:
    env.Prepend(CPPDEFINES=['SF_HAS_FFTW'])
```

and links `fftw` into the LIBS list for programs that need it. Follow the same pattern for other optional dependencies: guard with `env.get('LIBNAME')` and use `env.RSF_Place('sfprogname', None, var='LIBNAME', package='package-name')` to produce a stub binary when the library is absent.

---

## Choosing a language

Use this decision tree as a starting point:

- **C** — default choice for new contributions. The entire core library (`librsf`) is C, most existing programs are C, and contributions blend in naturally. Best for performance-critical code and programs that extend or wrap existing C library routines.
- **Python** — prototype quickly, leverage NumPy/SciPy for array math, or write research scripts that do not need the performance of a compiled binary. Python programs are installed as executable scripts (renamed from `M<name>.py` to `sf<name>`).
- **C++** — use when you need templates, object-oriented design, or integration with C++ libraries such as the Madagascar C++ API (`librsf++`) or Eigen/LAPACK wrappers.
- **Fortran 90** — preferred over F77 for new Fortran work. Use when the algorithm is naturally expressed in Fortran or when you are porting existing scientific Fortran code.
- **Fortran 77** — legacy only; avoid for new programs. Rarely used; no existing `M*.f` main programs in this tree — the per-language skill will point you at the `api/f77/Test*.*` files as templates.
- **Julia** — when your ecosystem is Julia or your collaborators prefer it; Julia programs are installed as scripts. Rarely used; no existing `M*.jl` main programs in this tree — the per-language skill will point you at the `api/julia/Test*.*` files as templates.
- **Matlab / Octave** — when the algorithm is already in `.m` form and conversion is not worth the effort; the build system uses configure flags to distinguish Matlab vs. Octave. Rarely used; no existing `M*.m` main programs in this tree — the per-language skill will point you at the `api/matlab/Test*.*` or `api/octave/Test*.*` files as templates.
- **Java** — rarely used; no existing `M*.java` main programs in this tree — the per-language skill will point you at the `api/java/Test*.*` files as templates.
- **Chapel** — parallel computing experiments; follow `user/*/M*.chpl` patterns.

---

## Next step

After loading this skill, load the language-specific skill for your target language:

```
skills/authoring-sf-programs-c/SKILL.md       # C
skills/authoring-sf-programs-python/SKILL.md  # Python
skills/authoring-sf-programs-cc/SKILL.md      # C++
skills/authoring-sf-programs-f90/SKILL.md     # Fortran 90
skills/authoring-sf-programs-julia/SKILL.md   # Julia
```

The language skill shows the exact boilerplate, `sf_init` / `sf_input` / `sf_output` call patterns, parameter retrieval idioms, and a worked example for that language.
