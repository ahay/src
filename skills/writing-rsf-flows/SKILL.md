---
name: writing-rsf-flows
description: Use when writing or modifying an SConstruct file that drives a Madagascar data-processing flow (Flow/Plot/Result/Fetch/Command).
---

## When to use

Use this skill whenever you are authoring or editing an `SConstruct` file that drives a Madagascar geophysical data-processing pipeline. This covers creating new processing workflows, adding stages to existing flows, configuring plots or results, fetching remote datasets, and integrating non-Madagascar steps via `Command`/`Action`. Any SConstruct that starts with `from rsf.proj import *` is using this DSL.

## The DSL surface

`framework/rsf/proj.py` exposes six top-level functions that wrap a singleton `Project` object (a subclass of SCons `Environment`). These are the only functions you should call directly in a user-facing SConstruct.

| Function | One-line description |
|----------|----------------------|
| `Flow(target, source, flow, **kw)` | Run a Madagascar command pipeline to produce one or more `.rsf` files. |
| `Plot(target, source, flow=None, **kw)` | Run a vplot command to produce a `.vpl` file (not a deliverable). |
| `Result(target, source, flow=None, **kw)` | Same as `Plot`, but copies the output into `Fig/` and registers it as a deliverable viewable with `scons view`. |
| `Fetch(file, dir, private=0, **kw)` | Download a data file from a remote server into the current directory. |
| `End(**kw)` | Finalize the build: wire up `view`, `print`, `lock`, and `test` aliases. Must be called at the end of every SConstruct. |
| `Command(targets, sources, action)` | SCons built-in re-exported; use for steps that are not Madagascar `sf*` programs. |

Full module-level signatures from `proj.py`:

```python
def Flow(target, source, flow, **kw): ...
def Plot(target, source, flow=None, **kw): ...
def Result(target, source, flow=None, **kw): ...
def Fetch(file, dir, private=0, **kw): ...
def End(**kw): ...
```

`Project.Flow` (the underlying method) has this full signature:

```python
def Flow(self, target, source, flow, stdout=1, stdin=1, rsfflow=1,
         suffix=sfsuffix, prefix=sfprefix, src_suffix=sfsuffix,
         split=[], np=1, reduce='cat', jobmult=1, local=0, noderotate=1,
         workdir=None, wall=''):
```

`Project.Fetch` full signature:

```python
def Fetch(self, files, dir, private=None, server=dataserver, top='data', usedatapath=True):
```

## Flow in depth

`Flow(target, source, flow)` is the workhorse of the DSL.

**Targets and sources map to `.rsf` files.** If `target` is the string `'spike'`, the output file is `spike.rsf`. If `source` is `'spike'`, the input is `spike.rsf`. The `.rsf` suffix is appended automatically (controlled by the `suffix` and `src_suffix` keyword arguments, both defaulting to `'.rsf'`).

**The `sf` prefix is implicit inside a `flow` string.** When `rsfflow=1` (the default), the build system prepends `sf` to each command token that is recognized as a Madagascar program. So you write `'spike n1=1000'` not `'sfspike n1=1000'`. This means you should never write the `sf` prefix inside a `Flow` command string.

**Pipes.** Multiple stages are separated with `|` inside the flow string, just like shell pipes:

```python
Flow('filtered_spike', None, 'spike n1=1000 k1=300 | bandpass fhi=2 phase=y')
```

**Ampersand-separated commands.** `&&` inside a Flow command string separates two sequential pipeline groups. Each group is assembled independently (with `sf` prefixes added, inputs/outputs wired) but they run in the same shell action in sequence, sharing the target output. Use it when you need two related commands that both contribute to producing the target — but if they're independent, prefer two separate `Flow` calls.

**Multi-target flows.** Pass a Python list as `target` to produce several output files from a single command. Use bare names (without `.rsf`) — the DSL appends the suffix automatically:

```python
Flow(['cmpt', 'offset'], 'input', 'some_program ...')
# or equivalently, a space-separated string:
Flow('out1 out2', 'input', 'program ...')
```

Note: names that already contain a `.` (e.g. `'cmpt.rsf'`) also work — the DSL skips suffix-adding when a `.` is already present — but bare names are the idiomatic form.

**Zero-input flows (`source=None`).** Pass `None` (not an empty string) when a command generates data with no RSF input:

```python
Flow('spike', None, 'spike n1=1000 k1=300')
```

Internally, when `source` is falsy the DSL sets `stdin=0`, meaning no RSF file is piped to stdin.

**Multiple sources.** Pass a space-separated string or a list. Inside the `flow` string, refer to extra sources as `${SOURCES[1]}`, `${SOURCES[2]}`, etc.:

```python
Flow('nmo', 'cmp offset vnmo',
     'nmo offset=${SOURCES[1]} velocity=${SOURCES[2]} half=n')
```

## Plot vs. Result

`Plot` and `Result` share the same signature:

```python
Plot(target, source, flow=None, **kw)
Result(target, source, flow=None, **kw)
```

**Two-argument shorthand.** When called with two positional arguments, the second is treated as `flow` and `source` is set to `target`:

```python
Plot('spike', 'wiggle clip=0.02 title="Spike"')
# equivalent to:
Plot('spike', 'spike', 'wiggle clip=0.02 title="Spike"')
```

**Plot** produces `<target>.vpl` in the current directory. It is an intermediate artifact — not tracked as a deliverable.

**Result** is identical to `Plot` but additionally:
- Copies the `.vpl` file into `Fig/<target>.vpl`.
- Registers the target for `scons view` (opens all results with `sfpen`).
- Registers the target for `scons print` and `scons lock`.
- Records the target name in `project.rest` for the `.rsfproj` manifest.

Rule of thumb: use `Plot` for intermediate composites or overlays; use `Result` for the final figure you want to present.

**Overlay example** (combining two plots):

```python
Plot('cmp', 'grey title="Synthetic CMP"')
Plot('time', 'graph yreverse=y wanttitle=n plotfat=3')
Result('cmp-time', 'cmp time', 'Overlay')
```

The second argument `'cmp time'` lists source `.vpl` files; `'Overlay'` is a built-in combine mode (see the `combine` dict in `proj.py`).

## Fetch

```python
Fetch(file, dir, private=0, **kw)
# Underlying Project.Fetch signature:
# Fetch(self, files, dir, private=None, server=dataserver, top='data', usedatapath=True)
```

`Fetch` downloads one or more files from a remote server. With the default `usedatapath=True`, the file is stored under `$DATAPATH` and a symlink is placed in the current directory (matching RSF's header-and-binary split). Key keyword arguments:

- `server` — URL of the server (default: `https://ahay.org` via `get_dataserver()`).
- `top` — top-level path component on the server (default: `'data'`).
- `private` — dict with `login`/`password`/`server` for FTP-authenticated downloads.
- `usedatapath` — if `True` (default), stores the binary in `DATAPATH` and symlinks.

Worked example from `book/rsf/tutorials/nmo/SConstruct`:

```python
Fetch('synthetic_cmp.npz', 'data',
      server='https://raw.githubusercontent.com',
      top='seg/tutorials-2017/master/1702_Step_by_step_NMO')
```

This downloads:
`https://raw.githubusercontent.com/seg/tutorials-2017/master/1702_Step_by_step_NMO/data/synthetic_cmp.npz`

After `Fetch`, the file is available as a source to subsequent `Flow` or `Command` calls.

## Command and Action

`Command` (re-exported from SCons) is for steps that are not `sf*` Madagascar programs — calling a Python function inline, running a third-party tool, or converting between file formats.

**Python function as an action:**

```python
def npz2rsf(target=None, source=None, env=None):
    import numpy, m8r
    data = numpy.load(str(source[0]))
    out = m8r.Output(str(target[0]))
    out.put('n1', data['arr'].shape[0])
    out.write(data['arr'])
    out.close()
    return 0

Command(['out.rsf'], 'input.npz', action=Action(npz2rsf))
```

The function signature `(target, source, env)` is required by SCons. `target` and `source` are lists of SCons `Node` objects; convert with `str(source[0])`.

**Shell command with a non-Madagascar tool:**

```python
Command('out.csv', 'in.rsf',
        'sfdd form=ascii < $SOURCE > $TARGET')
```

Use `$SOURCE` / `$TARGET` (SCons variables) for single-file targets inside a `Command` action string.

Do NOT use `Command` as a replacement for `Flow` when the command is an `sf*` program — use `Flow` so the DSL handles path resolution, binary file cleanup, and parallel execution correctly.

## Parameter interpolation

SConstruct files are plain Python, so you can compute parameters at build-configuration time and interpolate them into flow strings.

**`%` formatting (most common in existing Madagascar code):**

```python
fhi = 4
Flow('filter', 'spike', 'bandpass fhi=%g phase=y' % fhi)
```

**`.format()`:**

```python
Flow('filter', 'spike', 'bandpass fhi={fhi} phase=y'.format(fhi=fhi))
```

**f-strings (Python 3.6+):**

```python
Flow('filter', 'spike', f'bandpass fhi={fhi} phase=y')
```

Note: f-strings work fine in modern Madagascar (Python 3). The historic caution was that SCons uses `$VAR` substitution in command strings, so a literal `$` inside a flow string could clash. Because Madagascar flow strings are processed by `rsf.flow.Flow` before being handed to SCons, this is not usually a problem — but avoid `$` characters in f-string expressions just to be safe.

**Multiple parameters:**

```python
t0, v0 = 0.2, 4000
Flow('time', 'offset',
     'math output="sqrt(%g+input*input/%g)"' % (t0*t0, v0*v0))
```

## Running a flow

| Command | Effect |
|---------|--------|
| `scons` | Build all default targets (everything registered with `Flow`, `Plot`, `Result`). |
| `scons view` | Open all `Result` outputs with `sfpen` (requires a display). |
| `scons -n` | Dry run: print the commands that would be executed without running them. |
| `scons <target>` | Build only `<target>.rsf` (or `<target>.vpl`). Example: `scons filter`. |
| `scons print` | Print all results to the configured printer via `pspen`. |
| `scons lock` | Copy results to the locked figures directory for regression testing. |
| `scons -j4` | Parallel build using 4 jobs. |

The `scons view` alias is only registered if at least one `Result` call has been made.

## How it plugs into the larger build

Each chapter/tutorial lives in its own subdirectory with its own `SConstruct`. The top-level `book/SConstruct` and intermediate `book/<chapter>/SConstruct` files use `from rsf.book import *` instead of `from rsf.proj import *`. The `Book()` helper in `framework/rsf/book.py` recurses into subdirectories and compiles the LaTeX papers alongside the data flows.

A typical book hierarchy:

```
book/
  tutorial/
    SConstruct          # uses rsf.book — recurses into chapters
    users/
      SConstruct        # uses rsf.proj — data flows here
    authors/
      SConstruct
```

When writing data-processing flows (not book structure), always use `from rsf.proj import *`. Only use `from rsf.book import *` when you are building a multi-chapter book document. See `framework/rsf/book.py` for the `Book()` and `Paper()` helpers.

## Common mistakes

1. **Forgetting `End()`.**
   Every SConstruct that uses `rsf.proj` must call `End()` at the bottom. Without it, the `view`, `print`, `lock`, and `test` SCons aliases are never registered, and the `.rsfproj` manifest is never written.

2. **Writing `sfbandpass` instead of `bandpass` inside a Flow string.**
   The `sf` prefix is added automatically by the DSL — `Project.Flow` accepts a `prefix=sfprefix` parameter and prepends it to each command token. Writing `sfbandpass` explicitly is redundant and harder to read. It also reduces portability: if `sfprefix` were ever reconfigured (the knob exists, even if it's not the default), explicitly-prefixed commands would need to be updated everywhere. In practice `sfbandpass` still works because `add_prefix` in `framework/rsf/flow.py` checks whether the prefix is already present before adding it — so existing code using `sfbandpass` is not broken, just non-idiomatic. Stick with the bare form (`bandpass`) inside Flow strings.

3. **Using shell redirection instead of Flow pipelines.**
   Do not write `'spike n1=100 > out.rsf'`. Use `Flow('out', None, 'spike n1=100')`. Shell redirection bypasses the DSL's path resolution, binary-file tracking, and cleanup logic.

4. **Mixing relative and absolute paths in `source`.**
   Sources should be base names without the `.rsf` suffix (e.g., `'spike'`, not `'./spike.rsf'`). The DSL appends `.rsf` automatically. Passing full paths can confuse the dependency tracker.

5. **Calling `Result('name')` with only one argument.**
   `Result` requires at least two positional arguments: `Result(target, source_or_flow)`. With two arguments the second is treated as the `flow` string and source defaults to `target`. With three arguments: `Result(target, source, flow)`.

6. **Forgetting `source=None` for zero-input commands.**
   Commands like `sfspike` or `sfecho` read nothing from stdin. Pass `source=None` explicitly so the DSL sets `stdin=0`. Passing an empty string `''` is not equivalent.

7. **Using Python `print()` to inspect RSF headers at build time.**
   Use `sfin` or `sfattr` as separate shell commands, not inside the SConstruct. At SConstruct-read time the RSF files may not exist yet.

## Example

See `references/example-flow.SConstruct` for a self-contained runnable demo. Stage-by-stage walkthrough:

```python
from rsf.proj import *
```
Imports the DSL and initialises the singleton `Project`. This line is mandatory.

```python
Flow('spike', None, 'spike n1=1000 k1=300')
```
Runs `sfspike n1=1000 k1=300 > spike.rsf`. `source=None` because `sfspike` produces data from scratch.

```python
Flow('filter', 'spike', 'bandpass fhi=2 phase=y')
```
Runs `sfbandpass fhi=2 phase=y < spike.rsf > filter.rsf`. Note: `bandpass` not `sfbandpass`.

```python
fhi = 4
Flow('filter2', 'spike', 'bandpass fhi=%g phase=y' % fhi)
```
Demonstrates `%`-format interpolation. The value of `fhi` is baked in at SConstruct-parse time.

```python
Plot('filter', 'wiggle clip=0.02 title="Filtered spike"')
Plot('filter2', 'wiggle clip=0.02 title="Filter (fhi=%g)"' % fhi)
```
Each `Plot` produces a `.vpl` in the current directory. Two-argument form: source = target = `'filter'`.

```python
Result('filter', 'wiggle clip=0.02 title="Filtered spike"')
```
Same as the `Plot` above, but also writes `Fig/filter.vpl` and registers `filter` for `scons view`.

```python
End()
```
Wires up the `view`/`print`/`lock`/`test` aliases. Always the last line.

Running `scons` in the directory produces `spike.rsf`, `filter.rsf`, `filter2.rsf`, `filter.vpl`, `filter2.vpl`, and `Fig/filter.vpl`. Running `scons view` opens `Fig/filter.vpl` with `sfpen`.
