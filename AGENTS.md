# AGENTS.md — Madagascar

> Agent orientation for the Madagascar geophysical data analysis package.
> Authored for Claude Code, Copilot CLI, Codex, Gemini CLI. Also readable as plain project context by Cursor, Aider, Zed, and other agent harnesses.

## What Madagascar is

Madagascar is an open-source package for reproducible geophysical data analysis. Research and processing are expressed as SConstruct files that chain small composable `sf*` programs; each program reads and writes files in the RSF format (regularly sampled data, text header + binary).

## Mental model

1. **RSF files** are the data plane: one text header (`.rsf`) pointing at a binary (usually in `$DATAPATH`).
2. **`sf*` programs** are the compute plane: ~1,000 small Unix-style filters that read RSF on stdin and write RSF on stdout.
3. **SConstruct flows** are the orchestration plane: Python (`rsf.proj`) DSL using `Flow`, `Plot`, `Result`, `Fetch`, `Command`.
4. **vplot (`.vpl`)** is the visualization plane: `sfgrey`, `sfgraph`, `sfwiggle`, etc. produce `.vpl` files that `sfpen` displays.

## Repo map

- `api/` — language bindings. C is the reference; 9 other languages wrap it (C++, Python, F77, F90, Java, Julia, Matlab, Octave, Chapel).
- `framework/rsf/` — the Python DSL for flows: `proj.py` (Flow/Plot/Result/Fetch), `flow.py`, `tex.py`.
- `user/<author>/M<name>.[c|py|cc|f|f90|...]` — user-contributed main programs.
- `plot/main/` — vplot rendering programs (`grey.c`, `graph.c`, `wiggle.c`, `contour.c`, `dots.c`, `bargraph.c`, …).
- `book/` — 1,751 example SConstructs: teaching material, tutorials, published research reproductions.
- `system/` — filesystem-level helpers (`sfrm`, `sfcp`, …).

## RSF data model (briefly)

- Text header (`file.rsf`) contains `n1=..., d1=..., o1=..., label1=..., unit1=...` through axis 9.
- `n1` is the fastest-varying axis.
- Binary data is in a separate file pointed to by `in=` in the header, typically under `$DATAPATH` (set via `env.sh`).
- `<` redirects stdin from an RSF file; `>` writes to a new RSF (with its binary placed in `$DATAPATH`).

## Running things

- `scons` in any `book/` directory builds that directory's flow.
- `sfdoc <progname>` prints program help. Also `<progname>` with no args.
- `sfin file.rsf` inspects header contents.
- `sfattr < file.rsf` prints data statistics.
- `scons view` runs the flow and opens all `Result` plots.

## Discovering programs

- `sfdoc <name>` — full docs for a single program.
- `sfdoc -k <keyword>` — search docs.
- `user/<author>/M<name>.*` — read the source for how a program is built.

## Adding a program

Create a `M<name>.<ext>` file in `user/<youruser>/`. The `sf<name>` binary is produced automatically by the local `SConstruct`. See the `authoring-sf-programs` skill family for details per language.

## Skills

See `skills/` in this repo. Each path below resolves to a `SKILL.md` that the Claude Code (and compatible) Skill tool can load on demand.

| Skill | When to use |
| --- | --- |
| `skills/writing-rsf-flows/` | Writing or modifying an SConstruct that orchestrates RSF processing. |
| `skills/using-sf-programs/` | Composing a pipeline from existing `sf*` programs. |
| `skills/working-with-rsf-data/` | Inspecting, modifying, or debugging RSF headers and binaries. |
| `skills/plotting-with-vplot/` | Producing or debugging a vplot visualization. |
| `skills/authoring-sf-programs/` | Shared conventions for authoring a new `sf*` program (language-agnostic). |
| `skills/authoring-sf-programs-c/` | Authoring in C (reference implementation). |
| `skills/authoring-sf-programs-python/` | Authoring in Python. |
| `skills/authoring-sf-programs-cpp/` | Authoring in C++. |
| `skills/authoring-sf-programs-f90/` | Authoring in Fortran 90. |
| `skills/authoring-sf-programs-f77/` | Authoring in Fortran 77. |
| `skills/authoring-sf-programs-java/` | Authoring in Java. |
| `skills/authoring-sf-programs-julia/` | Authoring in Julia. |
| `skills/authoring-sf-programs-matlab/` | Authoring in Matlab. |
| `skills/authoring-sf-programs-octave/` | Authoring in GNU Octave. |
| `skills/authoring-sf-programs-chapel/` | Authoring in Chapel. |

## Gotchas

- `rm *.rsf` removes headers but leaves binaries behind in `$DATAPATH`. Use `sfrm` instead.
- Axis 1 is fastest-varying; `sftransp` swaps axes.
- Plot labels use escape codes: `\_` subscript, `\^` superscript, `\s<n>` size. Forgetting these in SConstruct strings leaves literal backslashes in output.
- The `sf` prefix is implicit inside a `Flow()` string — write `bandpass`, not `sfbandpass`.
- `$DATAPATH` must be writable and end in `/`.
- Complex data types (`sf_complex`) require `data_format=native_complex` (or `ascii_complex` / `xdr_complex`) in the header — the shorthand `complex` alone is not a valid value.
