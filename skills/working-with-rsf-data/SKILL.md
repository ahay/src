---
name: working-with-rsf-data
description: Use when inspecting, modifying, or debugging an RSF file, its headers, or its binary data placement.
---

## When to use

Use this skill for any task that touches RSF files directly:

- Inspecting what an RSF file contains (shape, type, axis metadata, value statistics).
- Modifying header fields (labels, units, axis parameters).
- Diagnosing "file not found", "wrong dimensions", or "garbled data" problems.
- Understanding why a pipeline produces unexpected results after copying or moving files.
- Cleaning up RSF files without orphaning their binaries.
- Working with complex-valued RSF data and choosing the right downstream program.
- Recovering from situations where headers and binaries have become separated.

## The two-file model

Every RSF dataset is split across two files:

1. **Text header** (`*.rsf`) — a plain-text key=value file you can read with `cat`. It describes axes, data type, and the path to the binary.
2. **Binary file** (`*.rsf@`) — the raw sample data, written under the directory named by `$DATAPATH`.

The `in=` field in the header connects the two. Here is a real header from `sfspike n1=100 n2=5 k1=50 mag=1 > data.rsf`:

```
4.3-git	sfspike	private/tmp/working-with-rsf-data-demo:	root@hostname	Sun Apr 19 12:28:16 2026

	d2=0.1
	o1=0
	n2=5
	o2=0
	label1="Time"
	data_format="native_float"
	label2="Distance"
	esize=4
	in="/var/tmp/data.rsf@"
	unit1="s"
	unit2="km"
	d1=0.004
	n1=100
```

Key observations:

- The first line is a provenance comment (program name, working directory, host, timestamp). It is not parsed by Madagascar programs.
- All key=value lines are indented with a leading tab character.
- `in=` holds the absolute path to the binary file. This path is set at creation time using `$DATAPATH`.
- `data_format="native_float"` is the default: 4-byte IEEE 754 float, native byte order.
- `esize=4` is the element size in bytes.

Because the binary lives under `$DATAPATH` and not next to the header, copying or moving `*.rsf` files with standard Unix tools (`cp`, `mv`, `rm`) only touches the header. The binary is untouched and can become orphaned.

## Axes

RSF supports up to nine axes, numbered 1 through 9. Each axis has five fields:

| Field     | Meaning                        | Default if omitted |
|-----------|--------------------------------|--------------------|
| `n<k>`    | number of samples on axis k    | (required for k=1) |
| `d<k>`    | sampling interval              | 1.0                |
| `o<k>`    | axis origin (value at index 0) | 0.0                |
| `label<k>`| axis label string              | "" (empty)         |
| `unit<k>` | unit string                    | "" (empty)         |

**Memory order**: axis 1 is the fastest-varying (innermost) axis in memory. For a 2D array `n1=100 n2=5`, sample `[i2][i1]` is at byte offset `(i2*100 + i1) * esize`. This is analogous to C's row-major layout where axis 1 is the column index.

In seismic usage, `n1` is commonly the time or depth axis (samples per trace) and `n2` and above are spatial or offset axes.

Trailing dimensions of size 1 are often omitted from `sfin` output. A file with `n1=100` and no `n2` is effectively 1-D.

When reading a header, if `d1` or `o1` are missing, treat them as 1.0 and 0.0 respectively. Programs that need axis metadata will use these defaults.

## Data types

The `data_format` field controls how the binary is encoded. The default is native 32-bit float.

**Float (default)**

```
data_format="native_float"
esize=4
```

This is what `sfspike`, `sfbandpass`, and most Madagascar programs produce. It is IEEE 754 single precision, native byte order of the machine that wrote the file. The `xdr_float` variant is big-endian XDR (used for cross-platform exchange). `ascii_float` writes human-readable ASCII values.

**Complex**

Complex-valued files use one of three prefixed format strings:

- `data_format="native_complex"` — pairs of native-float (re, im), 8 bytes per sample
- `data_format="ascii_complex"` — ASCII-encoded complex pairs
- `data_format="xdr_complex"` — XDR big-endian complex pairs

IMPORTANT: `data_format="complex"` is NOT a valid value. It is a common but incorrect shorthand. Always use one of the three prefixed forms above. Programs that check `data_format` will not recognize the bare string.

When `esize=8`, the file holds complex samples. Most float-consuming programs (`sfgrey`, `sfwiggle`, `sfattr`) do not accept complex input directly. Convert first:

- `sfcabs` — complex → float magnitude
- `sfreal` — complex → float real part
- `sfimag` — complex → float imaginary part

**Integer and byte**

```
data_format="native_int"    esize=4   # 32-bit signed integer
data_format="native_short"  esize=2   # 16-bit signed integer
data_format="native_uchar"  esize=1   # unsigned byte (type=byte in sfin output)
```

`sfin` reports `type=float`, `type=complex`, `type=int`, `type=short`, or `type=byte` based on `data_format`.

## Essential tools

**`sfin`** — displays the shape and type of one or more RSF files. Run it as `sfin file.rsf` (not via stdin). Output includes the `in=` path, `esize`, `type`, `form` (native/xdr/ascii), all axis parameters, and the total element and byte count. Use `sfin info=n file.rsf` to suppress axis detail and show only the binary path. This is the first tool to reach for when you encounter an unfamiliar file.

**`sfattr`** — summarizes the numeric content of an RSF file. Read from stdin: `sfattr < file.rsf`. Reports rms, mean, 2-norm, variance, standard deviation, max (with index), min (with index), nonzero sample count, and total sample count. Use `want=max` or `want=rms` to print a single statistic. Useful for sanity-checking that a pipeline produced reasonable values. Does not work on complex files without prior conversion.

**`sfput`** — updates header fields and writes a new output file (header + binary copy). Takes input on stdin and writes output on stdout: `< in.rsf sfput key1=val1 key2=val2 > out.rsf`. The output file gets its own independent binary. Use this to relabel axes, fix sampling intervals, or correct missing metadata. Any key that appears in an RSF header can be set with `sfput`. To remove a key, set it to the empty string. Note: because `sfput` copies the binary, it has disk and time cost proportional to file size.

**`sfrm`** — removes RSF files together with their binaries. Use it instead of Unix `rm`. Syntax mirrors `rm`: `sfrm file1.rsf file2.rsf`. Flags `-v` (verbose), `-f` (force, ignore missing), `-i` (interactive) are supported. Under the hood, `sfrm` reads the `in=` field from each header, deletes the binary at that path, then deletes the header. If the binary is already missing, `sfrm` still deletes the header. See also `sfmv` and `sfcp` for moving and copying RSF files correctly.

## Binary location

The `$DATAPATH` environment variable controls where Madagascar writes binary files. It must end with a trailing slash. When you run:

```bash
sfspike n1=1000 > spike.rsf
```

Madagascar writes the text header to `spike.rsf` in the current directory and the binary to `$DATAPATH/spike.rsf@`. The `in=` field in the header records the full absolute path.

**Why `rm *.rsf` leaks binaries**: deleting the text header with Unix `rm` removes only the key=value file. The binary under `$DATAPATH` has no back-reference to the header and becomes an orphan. On a busy workstation with `$DATAPATH=/var/tmp/`, orphaned binaries can accumulate silently.

Always use `sfrm` to delete RSF files. If you have already used `rm` and need to clean up orphans, list `$DATAPATH` and look for `*.rsf@` files that no longer have a matching header anywhere on the filesystem.

**Setting `$DATAPATH`**: the environment script `$RSFROOT/share/madagascar/etc/env.sh` sets a default `$DATAPATH`, often pointing to a system temp directory that may be periodically cleaned. If you need data to persist, set `DATAPATH` explicitly in your shell profile to a stable directory with sufficient space and a trailing slash, for example:

```bash
export DATAPATH=/home/user/rsf-data/
```

## Header manipulation

`sfput` writes a new header AND a new binary copy — it is NOT zero-cost for large files. The convenience is the compact syntax for updating header fields (labels, units, axis deltas). For metadata-only tweaks on multi-GB files, plan for the disk and time cost.

```bash
# Fix wrong sampling interval and add axis labels
< wrong.rsf sfput d1=0.002 label1="Depth" unit1="km" > corrected.rsf
```

`corrected.rsf` has its own independent binary; deleting `wrong.rsf` with `sfrm` does not affect it.

**`sfheadermath`** handles per-trace SEGY-style trace header fields. This is about per-trace header fields in segregated SEGY-style datasets, not about the `*.rsf` file's top-level metadata (which `sfput` manages). Not to be confused with the RSF file header. SEGY files have a binary header and per-trace headers stored in a separate `*.hdr` RSF file alongside the data `*.rsf`. `sfheadermath` reads those trace headers and computes new values from existing keys:

```bash
# Compute offset from shot and receiver x-coordinates
< data.hdr sfheadermath output="abs(sx-gx)" key="offset" > data_with_offset.hdr
```

Use `sfheadermath` when you need to derive or correct SEGY trace header quantities (offset, midpoint, elevation statics) without reprocessing the seismic samples themselves.

## Recovery scenarios

**Scenario A: missing binary**

Symptom: a program reports it cannot open or read the binary file.

1. Find the expected binary path: `grep "in=" file.rsf` or `sfin info=n file.rsf`.
2. Check whether the file exists at that path: `ls -lh /var/tmp/file.rsf@`.
3. If the file was created on a different machine or in an environment where `$DATAPATH` pointed to a now-inaccessible location, the `in=` path in the header will be unreachable. Use `sfput` to write a new file with `in=` pointing to the current, valid location — or copy the binary to the expected path: `< file.rsf sfput "in=/new/path/file.rsf@" > file_fixed.rsf`.
4. If the binary is truly gone (deleted with `rm` or temp-dir purge), the data is unrecoverable from the header alone.

**Scenario B: wrong axis labels or parameters**

Symptom: `sfin` shows incorrect `label1`, `unit1`, `d1`, or `o1`.

Use `sfput` to fix the labels. Because `sfput` writes an independent binary copy, `file_fixed.rsf` is safe to use on its own:

```bash
< file.rsf sfput label1="Time" unit1="s" d1=0.004 o1=0 > file_fixed.rsf
```

Verify with `sfin file_fixed.rsf`. Once confirmed correct, remove the original with `sfrm file.rsf` — this will not affect `file_fixed.rsf` since each has its own binary.

**Scenario C: segregated SEGY headers**

SEGY data imported with `sfsegyread` produces two files: `data.rsf` (sample data) and `data.hdr` (per-trace SEGY headers as a separate RSF file with integer type). Programs like `sfheadermath`, `sfheadersort`, and `sfintbin` consume the `.hdr` file. If the `.hdr` file is missing or mismatched, re-import with `sfsegyread` or reconstruct trace headers using `sfheadermath` from known geometry.

## Example

See `references/example-rsf-inspect.sh` for a self-contained walkthrough. Run it with:

```bash
bash references/example-rsf-inspect.sh
```

The script progresses through seven stages:

1. **Create a 2D RSF** — `sfspike n1=100 n2=5 k1=50 mag=1 > data.rsf` generates a 100-sample by 5-trace file with a spike at sample 50 in every trace.

2. **Cat the header** — `cat data.rsf` shows the plain-text key=value format. The `in=` line reveals the binary path under `$DATAPATH`.

3. **sfin** — `sfin data.rsf` prints the parsed shape: `n1=100 d1=0.004 o1=0 label1="Time"`, `n2=5 d2=0.1 o2=0 label2="Distance"`, element count, and byte count.

4. **sfattr** — `sfattr < data.rsf` reports rms=0.1, max=1 at position [50, 1], and that only 5 of 500 samples are nonzero (one spike per trace).

5. **sfput** — `< data.rsf sfput label2="Trace" unit2="" > labeled.rsf` writes a new header pointing to a new binary with updated metadata. `sfin labeled.rsf` confirms the change.

6. **Locate the binary** — `grep "^	in=" data.rsf` prints the `in=` line, showing the absolute path under `$DATAPATH`. This is what would be orphaned by `rm data.rsf`.

7. **sfrm** — `sfrm data.rsf labeled.rsf` removes both headers and both binaries. The final `ls` confirms no `*.rsf` files remain, demonstrating clean teardown.

Expected final output: `clean`.
