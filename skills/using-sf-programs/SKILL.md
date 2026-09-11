---
name: using-sf-programs
description: Use when composing a Madagascar data-processing pipeline from existing sf* programs — includes discovery, parameter conventions, and piping patterns.
---

## When to use

Use this skill any time you are assembling a Madagascar data-processing workflow
from existing `sf*` programs: creating synthetic data, filtering, reshaping,
transforming, inspecting, or plotting.  If you are *writing a new* `sf*` program
from scratch (in C or Python), consult the program-authoring skill instead.  If
you need to embed pipelines inside a reproducible SCons build, consult the
`writing-rsf-flows` skill.

The patterns here apply equally at the shell prompt, inside SConstruct `Flow`
calls, and inside short Python driver scripts — the programs themselves are
identical in all three contexts.

For a worked end-to-end pipeline that chains these programs to solve a named geophysical problem, see the `finding-workflow-examples` skill.

## The pipeline model

Every `sf*` program is a Unix filter.  It reads an RSF dataset from **stdin**
and writes an RSF dataset to **stdout**.  Because all programs speak the same
binary format, they compose without adapters:

```bash
# Inline pipe — no temporary files:
sfspike n1=500 n2=10 k1=200 mag=1 \
  | sfbandpass fhi=60 phase=n \
  | sfwindow n1=300 f1=100 \
  | sfattr
```

The equivalent using explicit intermediates:

```bash
sfspike n1=500 n2=10 k1=200 mag=1 > raw.rsf
< raw.rsf sfbandpass fhi=60 phase=n > filtered.rsf
< filtered.rsf sfwindow n1=300 f1=100 > windowed.rsf
< windowed.rsf sfattr
```

Both forms are semantically identical.  The pipe form is more concise for
throwaway work; the intermediate-file form is better when you want to
inspect or reuse mid-stream results.

Key properties of the model:

- **Composability is the design goal.** Each program does one thing and passes
  the header plus binary data downstream unchanged (except for the dimensions
  it modifies).
- **The header travels with the data.**  RSF files consist of a small ASCII
  header (`.rsf`) that records `n#`, `d#`, `o#`, `label#`, `unit#`, and a
  pointer to the binary data file.  Programs update the header automatically;
  you rarely touch it by hand.
- **Axis numbering is 1-based.**  `n1` is the fast (innermost) axis — usually
  time samples.  `n2` is the slow axis — usually trace index.  Higher axes
  follow naturally.
- **stdin/stdout are the default I/O names.**  Programs that take secondary
  inputs (e.g. a mask file) use named parameters: `mask=head.rsf`.

## Discovery

### Looking up a specific program

```bash
sfdoc sfbandpass
```

This prints the program description, synopsis with all parameters and defaults,
the known functions, and a list of reproducible examples that use the program.
Every parameter claim in this skill was verified with `sfdoc`.

### Keyword search

```bash
sfdoc -k bandpass
# → sfbandpass: Bandpass filtering.
# → sferf: Bandpass filtering using erf function.
# → sftrapepass: Trapezoid bandpass filter in the frequency domain.
```

`sfdoc -k <word>` searches program descriptions for the keyword and prints
matching program names and one-line summaries.  Use it when you know what you
want to do but not which program does it.

### Listing all sf* programs

```bash
ls "$RSFROOT/bin" | grep '^sf'
```

This lists every installed program.  Pipe through `grep` to narrow by name
fragment (e.g. `grep fft`).  The catalog is large (hundreds of programs);
rely on `sfdoc -k` for semantic search rather than scanning the full list.

Do **not** use `sfdoc -l <dir>` to list programs — that flag writes LaTeX
documentation to files and is not a listing tool.

### Reading source

Source lives under `user/<author>/M<name>.c` (user-contributed) or
`system/main/<name>.c` / `system/generic/M<name>.c` (core programs).
The `SOURCE` field in `sfdoc` output names the file.  Reading source is the
authoritative way to understand edge cases and defaults not covered by the
synopsis.

## Parameter conventions

Parameter conventions are uniform across all `sf*` programs.  The following
are grounded in actual `sfdoc` output.

### Axis parameters (`n#`, `d#`, `o#`, `label#`, `unit#`)

From `sfdoc sfspike`:

```
n#=          size of #-th axis
d#=[0.004,0.1,0.1,...]    sampling on #-th axis
o#=[0,0,...]              origin on #-th axis
label#=[Time,Distance,Distance,...]  label on #-th axis
unit#=[s,km,km,...]       unit on #-th axis
```

- `n1=1000` — 1000 samples along the first (fast) axis.
- `d1=0.004` — 4 ms sample interval (the default for axis 1).
- `o1=0` — axis starts at time 0.
- `label1="Time"` — string labels appear in plots and in `sfin` output.

The `#` is a literal digit: `n1`, `n2`, `n3`, …  Programs that accept axis
parameters accept as many as the data has.

### Boolean flags

From `sfdoc sfbandpass`:

```
phase=n [y/n]   y: minimum phase, n: zero phase
verb=n  [y/n]   verbosity flag
```

Booleans are always `y` or `n`.  Never use `true`/`false` or `1`/`0`.

### Comma-separated lists (multi-valued parameters)

From `sfdoc sfspike`:

```
k#=[0,...]   spike starting position  [nsp]
mag=         spike magnitudes  [nsp]
```

When `nsp=2`, you supply two values separated by commas:

```bash
sfspike n1=1000 n2=20 nsp=2 k1=300,700 mag=1,0.5
```

This places two spikes on axis 1 at samples 300 and 700 with magnitudes 1.0
and 0.5 respectively.  The `[nsp]` annotation in the sfdoc synopsis means the
parameter repeats once per spike.

**Note:** `sfspike` reads `nsp` BEFORE `k#`/`mag`, so you MUST set `nsp=N`
explicitly when passing N-element lists; it is not inferred from list length.
With the default `nsp=1`, only the first element of each list is used and only
one spike is placed.

### File-valued parameters

Some programs accept secondary files by name.  From `sfdoc sfheaderwindow`:

```
mask=    auxiliary input file name
```

Usage:

```bash
< data.rsf sfheaderwindow mask=selection.rsf > windowed.rsf
```

The named file is opened separately; stdin/stdout handle the primary data flow.

### The `output=` expression parameter

From `sfdoc sfmath`:

```
output=   Mathematical description of the output
```

This is a string containing a mathematical expression.  Named input files
become variable names:

```bash
sfmath x=file1.rsf y=file2.rsf output='sin((x+2*y)^power)' > out.rsf
sfmath < file1.rsf tau=file2.rsf output='exp(tau*input)' > out.rsf
```

When stdin is supplied, it is available as `input`.  When producing data from
scratch (no stdin), set `nostdin=y` and specify `n1=`, `d1=`, `o1=` etc.

### Window parameters (`f#`, `n#`, `j#`, `min#`, `max#`)

From `sfdoc sfwindow`:

```
f#=(0,...)        window start in #-th dimension
n#=(0,...)        window size in #-th dimension
j#=(1,...)        jump (decimation) in #-th dimension
min#=(o1,o2,...)  minimum in #-th dimension
max#=...          maximum in #-th dimension
```

You can mix coordinate-based (`min1=`, `max1=`) and sample-based (`f1=`, `n1=`)
windowing.  Unspecified parameters default to keeping the full extent.

## Core catalog

The table below covers the programs used most often, grouped by purpose.  Run
`sfdoc <name>` for full parameter lists.

### Synthesize

| Program   | What it does |
|-----------|-------------|
| `sfspike` | Generate spikes, boxes, planes, or constant arrays.  The primary tool for making synthetic data.  Key params: `n1`, `n2`, `k1`, `mag`, `nsp`. |
| `sfmath`  | Evaluate a mathematical expression to create or transform data.  Accepts named file inputs as variables.  Supports trig, log, exp, abs, erf, complex functions. |
| `sfnoise` | Add (or replace with) random noise.  Key params: `var`, `range`, `mean`, `type` (y=normal, n=uniform), `rep` (replace instead of add), `seed`. |

### Reshape

| Program     | What it does |
|-------------|-------------|
| `sfwindow`  | Extract a sub-volume along any axis.  Supports sample-based (`f#`, `n#`, `j#`) and coordinate-based (`min#`, `max#`) selection. |
| `sfpad`     | Pad with zeros.  Key params: `beg#` (prepend), `end#` (append), `n#` (set output length directly). |
| `sftransp`  | Transpose two axes.  Default: `plane=12` swaps axes 1 and 2.  Use `plane=13` for axes 1 and 3, etc. |
| `sfput`     | Set or override header parameters in place.  Useful for correcting `d1`, `o1`, `label1`, etc. without touching data. |
| `sfcat`     | Concatenate datasets along an axis.  Key param: `axis=` (default 3). |

### Filter

| Program       | What it does |
|---------------|-------------|
| `sfbandpass`  | Butterworth bandpass along axis 1.  Key params: `flo`, `fhi`, `nplo=6`, `nphi=6` (filter poles), `phase` (y=minimum-phase, n=zero-phase). |
| `sfsmooth`    | Triangle smoothing along any axis.  Key params: `rect1`, `rect2`, … (smoothing radius in samples on each axis).  `repeat=` applies the filter multiple times. |

### Transform

| Program   | What it does |
|-----------|-------------|
| `sffft1`  | FFT along axis 1 (time → frequency).  Key params: `inv=n` (forward), `sym=n`, `opt=y` (auto-pad to efficient length). |
| `sffft3`  | FFT along any extra axis (default `axis=2`).  Input and output are complex.  Key params: `axis`, `pad=2` (padding factor), `sym=n`. |
| `sfcabs`  | Convert complex RSF to float RSF by computing the complex magnitude.  Use after `sffft1`/`sffft3` before float-consuming programs such as `sfgrey`.  Use `sfreal` for the real part only. |

### Inspect

| Program          | What it does |
|------------------|-------------|
| `sfin`           | Print RSF header fields: `n1`, `d1`, `o1`, `label1`, … plus data-file path, element size, and a quick zero-check. |
| `sfattr`         | Print amplitude statistics: rms, mean, 2-norm, variance, std dev, max, min, nonzero sample count. |
| `sfheaderwindow` | Select traces whose header key satisfies a mask.  Key params: `mask=` (integer RSF file, nonzero = keep), `inv=n`. |
| `sfheadermath`   | Apply math to trace headers.  Key params: `key=` (header key to replace), `output=` (expression), `segy=y`. |

### Plot (brief — see `plotting-with-vplot` skill for full coverage)

| Program    | What it does |
|------------|-------------|
| `sfgrey`   | Raster (density / image) plot.  Writes a vplot byte stream. |
| `sfgraph`  | Line/graph plot.  Writes a vplot byte stream. |
| `sfwiggle` | Wiggle-trace plot with filled lobes.  Writes a vplot byte stream. |

These write vplot streams, not RSF.  Use `sfpen`, `xtpen`, or `pdfpen` to
render them, or chain into `pspen` for PDF output.  Full parameter coverage
is in the `plotting-with-vplot` skill.

## Piping patterns

### Pattern 1: Synthesize → filter → plot

The simplest end-to-end workflow: generate a synthetic impulse, apply a
filter, and plot the result.

```bash
sfspike n1=500 d1=0.004 o1=0 label1=Time unit1=s k1=250 \
  | sfbandpass fhi=60 phase=n \
  | sfwiggle title="Bandpass demo" | xtpen
```

Stage-by-stage:
1. `sfspike n1=500 d1=0.004 … k1=250` — creates a single 500-sample trace
   (axis 2 defaults to 1 trace) with a spike at sample 250.  Spike positioning
   is 1-based.  Omitting `k1` (or setting `k1=0`) fills the **entire** trace
   with a constant of magnitude `mag` — not what you usually want for an
   impulse test.  The axis metadata (`d1`, `label1`, `unit1`) flows
   downstream.
2. `sfbandpass fhi=60 phase=n` — zero-phase low-pass at 60 Hz.  `sfbandpass`
   reads `d1` from the header to convert Hz to normalized frequency; you do
   not need to repeat it.
3. `sfwiggle` — renders as a wiggle plot.  Axis labels and title come from
   header metadata.  Pipe to a viewer (`xtpen`, `pspen > out.ps`, etc.).

### Pattern 2: FK spectrum

Compute the 2D Fourier transform (time → frequency, space → wavenumber) and
display power.

```bash
sfspike n1=512 n2=64 d1=0.004 d2=25 k1=256 \
  | sfnoise var=0.01 \
  | sffft1 \
  | sffft3 axis=2 \
  | sfcabs \
  | sfgrey title="FK spectrum" | xtpen
```

Stage-by-stage:
1. `sfspike` — 512 × 64 array with a spike at time sample 256 on all traces.
2. `sfnoise var=0.01` — add Gaussian noise (variance 0.01) so the FK spectrum
   is not just a line.
3. `sffft1` — FFT along axis 1 (time).  Output is **complex**, shape changes
   to `(n1/2+1) × 64`.
4. `sffft3 axis=2` — FFT along axis 2 (offset/space).  Both axes are now in
   the frequency/wavenumber domain.  Output remains **complex**.
5. `sfcabs` — converts complex RSF to float RSF by computing the complex
   magnitude (|re + i·im|).  This is the required step before any
   float-consuming program like `sfgrey`.  Use `sfreal` instead if you only
   want the real part.  Do **not** use `sfmath output='abs(input)'` here:
   on complex input `sfmath abs` returns complex output, which causes
   `sfgrey` to error with "Need float input".
6. `sfgrey` — raster plot of the amplitude spectrum.

### Pattern 3: Transpose-filter-transpose

`sfbandpass` always filters along axis 1 (the fast axis).  When you want to
filter along axis 2 instead, transpose before filtering, filter, then
transpose back.

```bash
< data.rsf \
  sftransp \
  | sfbandpass fhi=30 phase=n \
  | sftransp \
  > filtered_axis2.rsf
```

Stage-by-stage:
1. `sftransp` — swaps axes 1 and 2.  What was the trace axis (axis 2) is now
   the fast axis (axis 1), which `sfbandpass` operates on.
2. `sfbandpass fhi=30 phase=n` — filters what was originally the trace axis.
3. `sftransp` — swaps back.  Axis order is restored to the original.

This transpose-operate-transpose idiom works for any axis-1-only operator
(FFT, smooth, etc.) when you need to operate on a different axis.

### Pattern 4: Multi-file arithmetic

`sfmath` accepts multiple named input files:

```bash
sfmath x=model.rsf obs=field.rsf output='obs-x' > residual.rsf
sfmath x=residual.rsf output='x^2' | sfattr want=mean
```

Each named parameter other than `output`, `type`, `datapath`, and `out` is
treated as a variable in the expression.  When stdin is also present it is
available as `input`.  When creating data from scratch, pass `nostdin=y` and
specify axis parameters explicitly.

## When NOT to pipe

Piping is not always the right choice.

**Cache expensive intermediates.**  If stage 2 of your pipeline takes 10 minutes
to compute (e.g., a migration or FWI gradient), write it to disk and reference
it explicitly.  In a SConstruct `Flow`, that is the natural structure anyway:
each `Flow` call is one stage with a named output.

```python
Flow('migrated', 'shot_data', 'sfkirchhoff ...')
Flow('result', 'migrated', 'sfagc | sfbandpass fhi=60')
```

**Debug a broken pipeline.**  Drop `sfin` or `sfattr` after each stage to check
that dimensions and values look right (see the Debugging section below).
Intermediates make that easy.

**Pipelines with branches.**  A pipeline is a linear chain; if two downstream
consumers need the same intermediate, write it to disk.  The shell `tee`
command can help but adds complexity.

**Very large datasets.**  Programs that must load an entire axis into RAM
(e.g., `sftransp` with large `plane=`) may need `memsize=` tuning.  Working
with named files makes it easier to profile memory at each stage.

## Debugging a broken pipeline

### Insert `sfin` or `sfattr` after each stage

Break the pipe at the suspect stage and write an intermediate, then inspect:

```bash
sfspike n1=1000 n2=20 k1=300 mag=1 > raw.rsf
sfin raw.rsf           # prints n1, d1, o1, label1, data-file path
< raw.rsf sfbandpass fhi=4 phase=y > filtered.rsf
sfin filtered.rsf      # verify dimensions after bandpass
sfattr < filtered.rsf  # check amplitude statistics
```

`sfin` prints header fields and a quick zero-check; `sfattr` prints rms,
mean, 2-norm, max, min, and variance.  Note that `sfin` terminates the
stream — it cannot be inserted transparently mid-pipe.

### Check `$DATAPATH`

RSF header files (`.rsf`) contain a pointer to a binary data file.  By
default, `DATAPATH` is `/var/tmp/` (or whatever your installation sets).  If
you see "cannot open data file" errors, check:

```bash
echo $DATAPATH
sfin myfile.rsf   # shows the actual data file path
```

Mismatches between `$DATAPATH` and the path recorded in the header happen
when files are moved without updating the header.  Use `sfput` to correct the
path, or re-run the pipeline with the correct `DATAPATH` set.

### Check axis metadata with `sfput`

If a downstream program complains about missing `d1` or wrong axis size, use
`sfput` to inject the correct values:

```bash
< data.rsf sfput d1=0.004 o1=0 label1=Time unit1=s > data_fixed.rsf
```

`sfput` passes all data through unchanged; it only rewrites header fields.

### Use `sfheaderwindow` for trace-by-trace header diagnostics

When working with SEGY-derived data that has trace headers, use
`sfheaderwindow` to select a subset of traces and `sfheadermath` to inspect
or correct header values:

```bash
# Select traces where offset < 2000 m (stored in a mask file):
sfheadermath < data.rsf output='offset<2000' key=mask > mask.rsf
sfheaderwindow < data.rsf mask=mask.rsf > near_offset.rsf
```

## Example

See `references/example-pipeline.sh` for a self-contained runnable example.
Run it with:

```bash
bash references/example-pipeline.sh
```

### Stage-by-stage walkthrough

**Stage 1 — Synthesize:**
`sfspike n1=1000 n2=20 nsp=2 k1=300,700 mag=1,0.5` creates a `1000 × 20`
dataset with **two** spikes at samples 300 and 700 (magnitudes 1.0 and 0.5).
`nsp=2` must be set explicitly — `sfspike` reads `nsp` before `k1`/`mag` and
does not infer it from list length; without it only the first spike would be
placed.  Defaults `d1=0.004`, `d2=0.1` come from `sfdoc sfspike`.

**Stage 2 — Bandpass:**
`sfbandpass fhi=4 phase=y` applies a 6-pole minimum-phase low-pass at 4 Hz.
`sfbandpass` reads `d1` from the header to convert Hz to normalized frequency;
you do not repeat `d1` on the command line.

**Stage 3 — Window:**
`sfwindow n1=500 f1=250` extracts samples 250–749 along axis 1, yielding
a `500 × 20` dataset.  Axis 2 is unchanged.

**Stage 4 — Transpose:**
`sftransp` (default `plane=12`) swaps axes 1 and 2, giving `20 × 500`.
All header fields (`n#`, `d#`, `o#`, `label#`) are updated automatically.

**Stage 5 — Summarize:**
`sfattr` prints rms, mean, 2-norm, variance, std dev, max, min, and sample
counts.  The `at 1 94` notation in the output gives the 2D sample location
of the maximum.

**Single-pipe form:** The script's final section chains all five programs with
`|`.  Output is identical — RSF programs are stateless filters; results depend
only on data and parameters, not on whether intermediates were written to disk.
