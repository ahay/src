---
name: authoring-sf-programs-octave
description: Use when authoring a Madagascar sf* main program in GNU Octave.
---

## When to use

Use this skill when you want to write a Madagascar `sf<name>` program using
GNU Octave — the open-source, MATLAB-compatible interpreter — and you do not
have a MATLAB license.

**Disambiguation: `.m` files are shared by Octave and MATLAB.**  A file named
`M<name>.m` could target either runtime. The difference is entirely in how
Madagascar is configured and how the helpers are installed:

- **Octave path** (this skill): pure `.m` helper files; no MEX compilation
  required; enabled with `API=octave` at configure time; installs the
  `rsf_*.m` helpers to Octave's function path.
- **MATLAB path**: C MEX extensions are compiled with `mex`; enabled with
  `API=matlab`; requires a MATLAB license. See
  `skills/authoring-sf-programs-matlab/SKILL.md`.

Typical users of this path: researchers who want MATLAB-like array syntax and
existing `.m` algorithms without a proprietary license.

Always load the shared conventions skill alongside this one:
[skills/authoring-sf-programs/](../authoring-sf-programs/SKILL.md)

---

## Skeleton

A minimal Octave `M<name>.m` script. No MEX compilation is needed — the Octave
API is implemented as pure `.m` helper files (`rsf_create.m`, `rsf_dim.m`,
`rsf_par.m`) that are placed on Octave's function path at install time.

```octave
% Mscale.m — sf<name> skeleton for GNU Octave
% One-line description of what this program does.
%
% Usage: sfscale < input.rsf scale=2.0 > output.rsf

% Read mandatory scalar parameter (no default -> error if absent)
[scale, st] = rsf_par('in.rsf', 'scale', []);
if st.err; error(st.msg); end
if isempty(scale); error('scale= required'); end

% Read optional integer parameter with default
[niter, st] = rsf_par('in.rsf', 'niter', 100);
if st.err; error(st.msg); end

% Query dimensions of input file
[dims, st] = rsf_dim('in.rsf');
if st.err; error(st.msg); end
n1 = dims(1);
n2 = dims(2);   % 1 if the dataset is 1-D

% Create output header by cloning the input header
st = rsf_create('out.rsf', 'in.rsf');
if st.err; error(st.msg); end

% Read binary data (rsf_read / rsf_write are NOT yet in api/octave;
% use sfdd or a shell pipe to convert to a plain binary float array)
fid_in  = fopen('<&0', 'rb');          % stdin when invoked as sf program
data    = fread(fid_in, n1*n2, 'float32');
fclose(fid_in);

% Process
data = data * scale;

% Write
fid_out = fopen('>&1', 'wb');          % stdout
fwrite(fid_out, data, 'float32');
fclose(fid_out);
```

> **Note on stdin/stdout**: When Madagascar invokes a program via a pipeline,
> `stdin` carries the raw RSF binary and `stdout` must produce the raw binary
> for the next stage. Octave scripts can use `fopen('<&0','rb')` and
> `fopen('>&1','wb')` for this purpose on POSIX systems.

---

## API cheat sheet

All helpers live in `api/octave/`. They call underlying `sf*` command-line
tools via `system()` and return a `stat` struct with fields `stat.err`
(logical) and `stat.msg` (string).

| Function | Signature | Purpose |
|---|---|---|
| `rsf_create` | `stat = rsf_create(out_filename, arg2)` | Write an RSF header. `arg2` is either an existing `.rsf` filename (copies that header) or a numeric vector of dimensions (creates a new header via `sfcreate`). |
| `rsf_dim` | `[dims, stat] = rsf_dim(in)` | Return a vector of dimensions for the RSF file `in`, with trailing length-1 dimensions stripped. Calls `sffiledims parform=n` internally. |
| `rsf_par` | `[par, stat] = rsf_par(file, name, default)` | Read scalar parameter `name` from header `file`. Returns `default` when the key is absent. Calls `sfget parform=n` internally. |

**Error handling pattern** (use consistently):

```octave
[val, st] = rsf_par('in.rsf', 'n1', 1);
if st.err
    error('rsf_par failed: %s', st.msg);
end
```

**`rsf_read` / `rsf_write`** are not present in `api/octave/` (only in
`api/matlab/` as MEX entry points). Read/write the binary payload directly
via Octave's `fread`/`fwrite` or by piping through `sfdd`.

---

## Build integration

`api/octave/SConstruct` is currently empty (placeholder); the `.m` helpers are
plain Octave function files that require no compilation step.

At configure time (`framework/configure.py`, `octave()` function):

1. `WhereIs('octave')` is run; if found, `env['OCTAVE']` is set.
2. `WhereIs('mkoctfile')` is checked; `env['MKOCTFILE']` is set if found.
   `mkoctfile` is the Octave function compiler (analogous to `mex`), but
   **it is not required for the pure `.m` API** — only needed if you later
   add compiled oct-files (`.oct`). The plain `rsf_*.m` helpers install
   without it.
3. Configure is triggered by passing `API=octave` (or including `octave` in
   the comma-separated `API` list) to `scons`.

**Installing the helpers**: the `rsf_*.m` files must be on Octave's function
path before your script can call them. Typical approaches:

```bash
# Option A: add to OCTAVE_PATH environment variable
export OCTAVE_PATH=/path/to/RSFROOT/lib:$OCTAVE_PATH

# Option B: addpath() inside your script (useful for SConstruct flows)
octave --eval "addpath('/path/to/RSFROOT/lib'); Mscale(...)"
```

**New user programs** — no MEX compile step is needed. Place `M<name>.m` in
`user/<youruser>/` and invoke it from a `Flow()` via the `OCTAVE` env
variable (similar to the pattern in
`book/rsf/school2025/plots/SConstruct`):

```python
# In user/<youruser>/SConstruct
octave = env.get('OCTAVE')
if octave:
    Flow('output', 'input',
         '%s --eval "addpath(...); Mscale(\'${SOURCE}\', \'${TARGET}\', ...); exit;"'
         % octave, stdin=0, stdout=-1)
```

Because the script is interpreted, the SConstruct does not need a compile
step — only an install/path setup.

---

## Pointers

Files in `api/octave/` with one-line descriptions:

| File | Description |
|---|---|
| `rsf_create.m` | Write an RSF header to disk — copies an existing header or creates one from a dimension vector by calling `sfcreate`. |
| `rsf_dim.m` | Return the dimension vector of an RSF file by calling `sffiledims parform=n`; strips trailing length-1 dimensions. |
| `rsf_par.m` | Read a named scalar parameter from an RSF header by calling `sfget parform=n`; returns a caller-supplied default when the key is absent. |
| `SConstruct` | Placeholder (empty); no compilation targets are needed for the pure `.m` Octave API. |

---

## Shared conventions

File naming, self-documentation comment format, parameter conventions, error
handling, testing, and build integration patterns that apply to every
`sf<name>` program regardless of language are documented in:

[skills/authoring-sf-programs/](../authoring-sf-programs/SKILL.md)

Key reminders from that skill relevant to Octave programs:

- Name your file `M<name>.m` inside `user/<youruser>/`.
- The installed binary is called `sf<name>` (build system drops `M`, prepends
  `sf`).
- The `.m` extension is shared with MATLAB; configure time (`API=octave` vs.
  `API=matlab`) determines which runtime is used and whether MEX compilation
  is required.
- Self-documentation: add a leading `%` comment block describing the program
  and its parameters so `sfdoc sf<name>` works correctly.
