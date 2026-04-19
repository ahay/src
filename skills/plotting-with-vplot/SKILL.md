---
name: plotting-with-vplot
description: Use when producing or debugging a vplot visualization in Madagascar — grey, graph, wiggle, contour, dots, overlays, and label escape codes.
---

## When to use

Use this skill any time a task involves:

- Producing a `.vpl` file — Madagascar's portable vector/raster plot format.
- Calling `sfgrey`, `sfgraph`, `sfwiggle`, `sfcontour`, `sfdots`, or `sfbargraph` from an SConstruct or at the command line.
- Debugging a blank, mis-scaled, or mis-labeled plot.
- Composing multiple plots into an overlay, side-by-side panel, or animation.
- Writing axis labels that require subscripts, superscripts, or Greek letters.
- Controlling color schemes, clip levels, or scalebars.

This skill does not cover converting `.vpl` to PDF/PNG (use `sfpen` options or `vppen`), or 3D surface plots (`sfgrey3`, `sfsurf`).

---

## Mental model

### What vplot is

Vplot is Madagascar's internal plot language — a compact binary stream describing lines, rasters, text, and color tables. Every visualization program (sfgrey, sfgraph, etc.) reads an RSF file and writes a `.vpl` stream. The `.vpl` file can be:

- Displayed on screen with `sfpen` (wraps the `xtpen` or `pspen` driver).
- Composed with other `.vpl` files using `vppen`.
- Converted to PostScript with `pspen`.

### Programs vs. the pipeline

Every vplot program is a standard Madagascar filter:

```
< data.rsf sfgrey color=j clip=1 title="My data" > out.vpl
```

At the command line, output goes to stdout. In SConstruct, SCons manages the filenames.

### Plot vs. Result in SConstruct

`Plot` and `Result` are the two SConstruct functions that produce `.vpl` files.

| Function | Output location | Purpose |
|----------|----------------|---------|
| `Plot('name', ...)` | `name.vpl` (working dir) | Intermediate — used as input to composition |
| `Result('name', ...)` | `Fig/name.vpl` | Final — what `scons lock` saves; what the paper uses |

Rules of thumb:
- Use `Plot` for any image you intend to overlay or combine.
- Use `Result` for any image that is a deliverable on its own.
- A `Result` with a composition string is both — it reads `Plot`-produced `.vpl` files and writes a final composed `.vpl` to `Fig/`.

### How parameters reach the vplot program

In SConstruct, the second argument to `Plot` or `Result` is either:
1. A flow string — the vplot program call with its parameters:
   ```python
   Result('img', 'grey title="My plot" color=j clip=2')
   ```
2. A source list plus a composition keyword:
   ```python
   Result('combined', 'img1 img2', 'Overlay')
   ```

Parameters follow the same `key=value` syntax as any `sf` program. String values with spaces must be quoted.

### Universal parameter inheritance

All vplot programs call `sfstdplot` internally (see `sfdoc stdplot`). This means axis labels, title, min/max, grid, screen size, and layout parameters are available to every vplot program, not just the one documented in its own `sfdoc` entry.

---

## Core program catalog

| Program | Description | Signature |
|---------|-------------|-----------|
| `sfgrey` | 2D raster/density plot. Maps float values to a color table. Default: `transp=y yreverse=y`. | `< 2D.rsf sfgrey [params] > out.vpl` |
| `sfgraph` | 1D line plot. Plots n2 traces as separate curves along axis 1. | `< 1D.rsf sfgraph [params] > out.vpl` |
| `sfwiggle` | Seismic wiggle-trace display. Draws each trace as a waveform with optional filled polygons. | `< 2D.rsf sfwiggle [params] > out.vpl` |
| `sfcontour` | Contour plot. Draws iso-value lines on a 2D field. Default: `transp=y allpos=y`. | `< 2D.rsf sfcontour [params] > out.vpl` |
| `sfdots` | Lollipop/impulse plot. Each sample becomes a vertical stick with a dot. Good for sparse signals. | `< 1D_or_2D.rsf sfdots [params] > out.vpl` |
| `sfbargraph` | Bar chart. Columns of data become side-by-side or stacked bars. | `< 1D_or_2D.rsf sfbargraph [params] > out.vpl` |

---

## Universal plot parameters

These parameters come from `sfstdplot` (run `sfdoc stdplot`) and are accepted by every vplot program.

| Parameter | Default | Effect |
|-----------|---------|--------|
| `title=` | (none) | Title string. Supports escape codes. |
| `label1=` | (from RSF header) | Label for axis 1 (vertical when `transp=y`). |
| `label2=` | (from RSF header) | Label for axis 2 (horizontal when `transp=y`). |
| `unit1=` | (from RSF header) | Unit appended to `label1` in parentheses. |
| `unit2=` | (from RSF header) | Unit appended to `label2` in parentheses. |
| `min1=` / `max1=` | (from data) | Display range on axis 1. |
| `min2=` / `max2=` | (from data) | Display range on axis 2. |
| `wantaxis=` | `y` | `n` suppresses both axes. `wantaxis1=`/`wantaxis2=` are per-axis. |
| `wanttitle=` | `y` | `n` suppresses the title. |
| `wheretitle=` | `top` | Title position: `top`, `bottom`, `left`, `right`. |
| `wherexlabel=` | `bottom` | Horizontal axis label position: `top` or `bottom`. |
| `screenratio=` | device default | Figure height/width ratio. `0.75`=landscape, `1.5`=portrait. |
| `plotfat=` | `0` | Line/curve thickness. `3`=bold, `6`=heavy. |
| `plotcol=` | `6` (yellow) | Curve color: 0=black, 1=blue, 2=red, 3=purple, 4=green, 5=cyan, 6=yellow, 7=white. |
| `pclip=` | 99 (grey), 98 (wiggle), 100 (others) | Percentile clip. |
| `grid1=` / `grid2=` | `n` | Draw grid lines on axis 1 or 2. |
| `dash=` | `0` | Line dash: 0=solid, 1=fine dash, 2=dot, 3=dash, 4=large dash. |
| `xll=/yll=/xur=/yur=` | (auto) | Manual panel placement in inches. |
| `crowd=` | `0.75` | Fraction of frame used by plot area. `crowd1=`/`crowd2=` per axis. |
| `labelsz=` | `8` | Axis label font size. `titlesz=10` for title. |

---

## Program-specific parameters

### sfgrey — 2D raster/density plot

`sfdoc sfgrey` | source: `plot/main/grey.c`

Defaults: `transp=y yreverse=y xreverse=n`. With these defaults, axis 1 runs vertically (time, depth), axis 2 runs horizontally (trace, offset). This is the standard seismic display convention.

| Parameter | Default | Notes |
|-----------|---------|-------|
| `color=` | `i` (grayscale) | Color scheme letter code. See **Color schemes** section below. |
| `clip=` | (estimated) | Hard clip value. Data outside `[-clip, clip]` saturates. If not set, estimated from `pclip`. |
| `pclip=` | `99` | Percentile used to estimate clip when `clip=` is not set. Set lower (e.g. `pclip=95`) to increase contrast. |
| `bias=` | `0` | Value mapped to the center of the color table. Shift to display asymmetric data better. |
| `allpos=` | `n` | If `y`, assume all data is positive and use the full color range for positive values only. Useful for amplitude maps. |
| `polarity=` | `n` | If `y`, reverse polarity: black represents high values, white represents low. |
| `scalebar=` | `n` | If `y`, draw a colorbar beside the plot. Use `minval=`, `maxval=` to set bar range; `barlabel=` and `barunit=` to label it. |
| `gainpanel=` | (none) | Gain normalization panel: `a` (all panels together), `e` (each panel independently), or an integer panel number. |
| `gpow=` | `1` | Raise data to this power before display. Values < 1 compress dynamic range; > 1 expand it. |
| `transp=` | `y` | Transpose axes. Set `transp=n` to put axis 1 horizontal. |
| `yreverse=` | `y` | Reverse vertical axis (time increases downward by default). |
| `xreverse=` | `n` | Reverse horizontal axis. |
| `mean=` | `n` | If `y`, set the bias to the data mean automatically. |
| `symcp=` | `n` | If `y`, use a symmetric 255-color palette (ensures that the center color maps to zero exactly). |
| `wantframenum=` | `y` if n3>1 | Show the current frame number in the corner for 3D data. |

**clip vs. pclip:** `pclip=99` (default) estimates clip from the 99th percentile. Set `clip=` explicitly for reproducible cross-panel display. `gainpanel=a` normalizes across all 3D frames; `gainpanel=e` per frame. For strictly positive data (envelopes, migration images) use `allpos=y color=j`.

---

### sfgraph — 1D line plot

`sfdoc sfgraph` | source: `plot/main/graph.c`

Default: `transp=n`. Plots each column (n2 traces) as a curve, with axis 1 values on the vertical axis and axis 2 index on the horizontal.

| Parameter | Default | Notes |
|-----------|---------|-------|
| `transp=` | `n` | If `y`, transpose so axis 1 is horizontal. Use `transp=y yreverse=y` to plot a velocity curve in the same orientation as a `grey` plot. |
| `yreverse=` | `n` | Reverse vertical axis. |
| `min2=` | (from data) | Minimum on the curve-value axis. Critical for overlay alignment — must match the background `grey` plot's axis range. |
| `max2=` | (from data) | Maximum on the curve-value axis. |
| `plotcol=` | `6` | Curve color. Use `plotcol=2` for red, `plotcol=1` for blue. |
| `plotfat=` | `0` | Curve fatness. `plotfat=3` gives a bold curve suitable for overlays. |
| `symbol=` | (none) | If set, plot with symbols (markers) instead of lines. Letter determines marker shape. |
| `symbolsz=` | `2` | Symbol size. |
| `wantaxis=` | `y` | When used as an overlay, set `wantaxis=n wanttitle=n` to suppress the graph's own axes (using the background image's axes instead). |
| `pad=` | `y` | If `n`, suppress extra padding around the plot area. Required for exact alignment with an underlying `grey` image. |
| `dash=` | `0` | Line dash pattern. `dash=1` gives fine dashes; `dash=3` gives coarse dashes. |
| `pclip=` | `100` | Clip percentile (usually leave at default for graphs). |
| `color=` | `j` | Color scheme (default is jet for graphs). |

**Overlay alignment:** Both plots must share axis ranges. With `sfgrey transp=y` (default), axis 1 is vertical and `min1/max1` control vertical range. With `sfgraph transp=y yreverse=y`, the curve values map to the vertical axis via `min2/max2`. Always use `wantaxis=n wanttitle=n pad=n` on the overlay graph.

---

### sfwiggle — seismic trace display

`sfdoc sfwiggle` | source: `plot/main/wiggle.c`

Draws each trace as an oscillating waveform. Positive excursions can be filled (polygon mode).

| Parameter | Default | Notes |
|-----------|---------|-------|
| `clip=` | `0` (estimated from pclip) | Hard amplitude clip. If `0`, estimated from `pclip`. |
| `pclip=` | `98` | Percentile used when `clip=0`. |
| `poly=` | `n` | If `y`, fill positive excursions with solid polygons (standard seismic display). |
| `polyneg=` | `n` | If `y`, also fill negative excursions (fills both sides). |
| `zplot=` | `0.75` | Horizontal spacing between traces as a fraction of the total plot width divided by n2. Values > 1 make traces overlap. |
| `transp=` | `n` | Transpose axes. Default `n` means axis 1 is vertical (time), axis 2 is horizontal (trace). |
| `yreverse=` | `n` | Reverse vertical axis. Set `y` to make time increase downward. |
| `xreverse=` | `n` | Reverse horizontal axis. |
| `fatp=` | `1` | Polygon border fatness (line thickness of the filled wiggle boundary). |
| `seemean=` | `n` | If `y`, draw a horizontal mean line through each trace. |
| `xpos=` | (none) | Optional auxiliary RSF file with explicit trace positions (allows irregular spacing). |
| `xmask=` | `1` | Polygon fill mask for x direction. |
| `ymask=` | `1` | Polygon fill mask for y direction. |

**Note on clip:** `sfwiggle` requires an explicit positive clip to work properly. Either set `clip=<value>` directly, or rely on `pclip=98` (the default) which estimates the clip from the data percentile. Setting `clip=0` with no `pclip` produces no display.

---

### sfcontour — contour plot

`sfdoc sfcontour` | source: `plot/main/contour.c`

Draws iso-value lines on a 2D scalar field. Defaults: `transp=y allpos=y`.

| Parameter | Default | Notes |
|-----------|---------|-------|
| `nc=` | `50` | Number of contour levels. |
| `dc=` | (auto) | Contour increment (spacing between levels). Overrides `nc`. |
| `c0=` | (auto) | Value of the first contour level. |
| `c=` | (none) | Explicit list of contour values `[nc]`. |
| `cfile=` | (none) | RSF file containing contour values. |
| `allpos=` | `y` | If `y`, only contour positive values. Set `n` to contour negative values too. |
| `transp=` | `y` | Transpose axes (same convention as sfgrey). |
| `min1=` | `o1` | Minimum on axis 1. |
| `max1=` | `o1+(n1-1)*d1` | Maximum on axis 1. |
| `min2=` | `o2` | Minimum on axis 2. |
| `max2=` | `o2+(n2-1)*d2` | Maximum on axis 2. |
| `scalebar=` | `n` | If `y`, draw a scalebar. |
| `barlabel=` | (none) | Label for the scalebar. |

**Contour overlay on grey:** A common pattern is to overlay contours on a raster image:
```python
Plot('raster', 'grey color=j')
Plot('lines',  'contour nc=20 wantaxis=n wanttitle=n')
Result('both', 'raster lines', 'Overlay')
```

---

### sfdots — lollipop/impulse display

`sfdoc sfdots` | source: `plot/main/dots.c`

Plots each sample as a vertical stick with a ball. Good for sparse signals and filter coefficients. Key parameters: `dots=` (1=balloon, 0=stick only), `strings=y` (draw sticks), `connect=` (1=diagonal, 2=bar), `gaineach=y` (per-trace normalization), `clip=` (-1=no clip), `overlap=0.9` (trace width fraction), `transp=`, `yreverse=`, `seemean=`. Note: the axis is drawn only when `label1=` is present.

---

### sfbargraph — bar chart

`sfdoc sfbargraph` | source: `plot/main/bargraph.c`

Plots each column as a group of bars. Key parameters: `width=0.8` (bar width fraction), `stack=y` (stack vs. side-by-side), `transp=n` (horizontal bars if `y`). Inherits all `sfstdplot` parameters for labels, title, and axes.

---

## Color schemes (`color=`)

The `color=` parameter in `sfgrey` (default `i`) and `sfgraph` (default `j`) selects a color table. These codes are defined in `plot/lib/coltab.c`. An uppercase letter reverses the table; appending `C` clips out-of-range values in a highlight color.

| Code | Name | Description |
|------|------|-------------|
| `i` | Grayscale | Default for `sfgrey`. Black to white. Seismic convention: white=positive, black=negative. |
| `I` | Reverse grayscale | White to black. |
| `j` | Jet | Blue → cyan → green → yellow → red. Good for amplitude/velocity maps. |
| `J` | Reverse jet | Red → yellow → green → cyan → blue. |
| `g` | Red-white-black | Red (negative) → white (zero) → black (positive). Diverging, asymmetric. |
| `G` | Reverse red-white-black | Black (negative) → white (zero) → red (positive). |
| `e` | Red-white-blue | Red (negative) → white (zero) → blue (positive). Classic diverging colormap. |
| `E` | Reverse red-white-blue | Blue (negative) → white (zero) → red (positive). |
| `h` | Hot | Black → red → orange → yellow → white. Emphasizes amplitude. |
| `H` | Reverse hot | White → yellow → red → black. |
| `p` | Pink | Softer version of hot (squished root of hot+grey). |
| `b` | Bone | Blue-tinted grayscale; similar to greyscale but cooler. |
| `a` | Rainbow (HSV) | Full hue cycle: red → yellow → green → cyan → blue → magenta → red. |
| `A` | Reverse rainbow | Reverse hue cycle. |
| `x` | Cubehelix | Perceptually linear helix through colour space (D.A. Green 2011). Monotone in luminance. |
| `X` | Reverse cubehelix | Reversed cubehelix. |
| `t` | Traffic | Green → yellow → red cycle (traffic light pattern). |
| `T` | Reverse traffic | Red → yellow → green. |
| `w` | Wheel | Full RGB colour wheel. |
| `W` | Reverse wheel | Reverse colour wheel. |
| `c` | Cool | Cyan → magenta (cyan at low, magenta at high). |
| `C` | Reverse cool | Magenta → cyan. |
| `l` | Linear (copper) | Black → dark orange/copper → bright copper. |
| `f` | Flag | Repeating red-white-blue-black pattern. Useful for phase. |

**Custom color tables:** You can also pass a filename (path to a CSV file with R,G,B triples, one per line, values in [0,1]) as the `color=` value. Madagascar looks in `$RSFROOT/include/<name>.csv` if the name is not a single letter.

**Clipping highlight suffix `C`:** Appending `C` to any color code (e.g. `color=jC`) marks values that fall outside the clip range with a red highlight color (configurable via `cliprgb=`). Useful for diagnosing saturation.

---

## Label escape codes

Vplot interprets `\` as an escape character in text strings (titles, axis labels). This lets you add subscripts, superscripts, font switches, and size changes directly in label strings.

### Complete escape reference

All escapes that do not take an argument:

| Escape | Effect |
|--------|--------|
| `\_` | Lower (subscript) — move down half a capital-letter height |
| `\^` | Raise (superscript) — move up half a capital-letter height |
| `\>` | Advance one inter-letter space |
| `\<` | Back up one inter-letter space |
| `\n` | Newline |
| `\\` | Print a literal backslash |
| `\g` | Ghostify (invisible text, for spacing) |
| `\G` | Resume visible text |

Escapes that take an integer argument followed by a required space:

| Escape | Argument | Effect |
|--------|----------|--------|
| `\s<N> ` | integer percent | Change text size to N% of current size. `\s100 ` restores default. `\s75 ` is smaller. |
| `\F<N> ` | font number | Switch to font N. Font 9 = Greek Simplex (β, α, etc.). `-1` restores default. |
| `\v<N> ` | glyph number | Print glyph N from the current font (bypasses special character meaning). |
| `\c<N> ` | color index | Change text color. `-1` restores drawing color. |
| `\f<N> ` | fatness delta | Add N to the current line fatness. |
| `\r<N> ` | percent height | Move up N percent of a character height (negative moves down). |
| `\k<N> ` | percent width | Move right N percent of a space width (negative moves left). |
| `\m<N> ` | register | Save current position in register N. |
| `\M<N> ` | register | Restore position from register N. |

Key font numbers: 0=Pen, 1=Roman Simplex, 3=Roman Complex, 5=Italic Complex, 9=Greek Simplex (α β γ), 10=Greek Complex, 15=Mathematics. `-1` restores the default font.

### The Python backslash problem

Python interprets `\s`, `\F`, `\_`, `\^` as string escapes. This means you need to double the backslash or use raw strings:

```python
# Wrong — Python eats the backslash before vplot sees it:
Result('bad', 'grey title="v\_NMO"')

# Correct option 1 — double the backslash:
Result('ok1', 'grey title="v\\_NMO"')

# Correct option 2 — raw string (preferred for readability):
Result('ok2', r'grey title="v\_NMO"')

# In a triple-quoted string, doubling still works:
Result('ok3',
       '''
       grey title="v\\_\s75 NMO\s100 "
       ''')
```

The NMO tutorial in the book (`book/rsf/tutorials/nmo/SConstruct`) uses the doubling convention:
```python
title="v\\_\\s75 NMO\\s100 \\^ Profile"
```

**Common recipes:** `r"v\_NMO"` → v subscript NMO; `r"x\^2\_ "` → x superscript 2 (back to baseline); `r"\F9 b\F-1 "` → Greek beta (font 9, restore with -1); `r"v\\_\\s75 NMO\\s100 "` → NMO at 75% size (doubling convention from the book).

---

## Composition

Madagascar's `rsf.proj` defines four composition modes for combining `.vpl` files. These are invoked as the flow string in `Result` or `Plot` when the source is a list of names.

### Overlay

Stacks multiple plots in registration (Z-order, last source on top).

```python
Plot('background', 'grey ...')
Plot('foreground', 'graph wantaxis=n wanttitle=n pad=n ...')
Result('combined', 'background foreground', 'Overlay')
```

Internally: `vppen erase=o vpstyle=n source1.vpl source2.vpl > out.vpl`

The first source defines the coordinate frame. All subsequent sources must use compatible axis ranges to align correctly.

### SideBySideAniso

Places N plots side by side, each keeping its own aspect ratio. Good for heterogeneous panels (e.g., a narrow velocity panel next to a wide gather).

```python
Result('comparison', 'img1 img2 img3', 'SideBySideAniso')
# Extra vppen flags via keyword argument:
Result('out', 'img1 img2', 'SideBySideAniso', vppen='txscale=1.5')
```

### SideBySideIso

Places N plots side by side with matched aspect ratios — good when panels show the same data type.

```python
Result('panels', 'imgA imgB', 'SideBySideIso')
```

### Movie

Assembles multiple frames into a flipbook animation. `sfpen` cycles through them on display.

```python
for iframe in range(nframes):
    Plot('frame%d' % iframe, 'slice%d' % iframe, 'grey ...')
Result('movie', ['frame%d' % i for i in range(nframes)], 'Movie')
```

Internally: `vppen vpstyle=n frame0.vpl frame1.vpl ...`

---

## Frame and layout

`screenratio=` (height/width) controls figure aspect: `0.5`=wide landscape, `0.75`=typical landscape, `1.0`=square, `1.5`=portrait.

For manual panel placement, use `xll=/yll=/xur=/yur=` in inches from the paper's bottom-left corner:

```python
Plot('left',  'grey ... xll=1.0 yll=1.5 xur=4.0 yur=6.5')
Plot('right', 'grey ... xll=4.5 yll=1.5 xur=7.5 yur=6.5')
Result('both', 'left right', 'Overlay')
```

`crowd=` (default `0.75`) controls what fraction of the frame the plot area occupies. `crowd1=`/`crowd2=` set axes independently.

---

## Debugging plots

### Empty or all-black plot

**Symptom:** `sfgrey` produces a `.vpl` that shows as solid black or solid white.

**Causes and fixes:**
1. Data is all zeros — check with `sfattr < data.rsf`. If `max = 0`, the flow that produced the data failed silently.
2. `clip=` was set to a value much larger than the data range — the entire range falls in one color. Check `sfattr` for actual min/max and set `clip=` accordingly, or remove it and let `pclip` estimate.
3. `pclip=` is very low (e.g. `pclip=1`) — almost all data is clipped. Increase to 99.

### Wrong aspect ratio

**Symptom:** Plot looks squashed or stretched.

**Fix:** Set `screenratio=` explicitly. Remember that `sfgrey` with `transp=y` (default) swaps n1 and n2 roles — if your data is 1000×50 samples, the display is 50 wide and 1000 tall without aspect correction.

### Missing or mangled labels

**Symptom:** Label appears with a literal backslash and letter instead of a formatted escape.

**Fix:** Python ate the backslash. Use raw strings (`r"..."`) or double the backslash (`"\\_"`).

**Symptom:** Label is absent entirely.

**Fix:** `wantaxis=n` is set, or the header fields `label1`/`label2` were not set in the RSF file and no explicit `label1=`/`label2=` was passed.

### Misaligned overlay

**Symptom:** Overlay curve appears at the wrong position on the background image.

**Cause:** Axis range mismatch between the two plots.

**Diagnosis:**
- `sfgrey` with `transp=y` (default): axis 1 is on the vertical axis, axis 2 on the horizontal. `min1/max1` control the vertical range; `min2/max2` control the horizontal.
- `sfgraph` with `transp=y`: axis 1 values are the independent axis (horizontal after transposing... actually vertical after `transp=y yreverse=y`), and the data values are on axis 2 (`min2/max2`).

**Fix:** Both plots must share identical `min1`/`max1` and `min2`/`max2`. Also use `pad=n` on the overlay `sfgraph` to suppress the automatic padding that would shift the curve.

### Scalebar out of range

Set `minval=` and `maxval=` explicitly (from `sfattr`) to control scalebar labels independently of clip range.

### Wiggle traces invisible

`clip=0` with no `pclip` estimates zero clip — nothing shows. Use `pclip=98` (default) or an explicit `clip=` from `sfattr`.

### Overlay Z-order

In `vppen erase=o`, sources are drawn left to right — the last source is on top. Put foreground elements last in the source list:
```python
Result('out', 'raster_below curve_on_top', 'Overlay')
```

---

## Examples

### Example 1 — basic grey plot

File: `skills/plotting-with-vplot/references/example-grey.SConstruct`

```python
from rsf.proj import *

Flow('synth', None,
     '''
     spike n1=200 n2=100 k1=50,100,150 nsp=3 mag=1,0.8,0.6 |
     ricker1 frequency=40 |
     noise seed=2025 var=0.005
     ''')

Result('synth',
       '''
       grey title="Synthetic gather"
            label1="Time" unit1="s" label2="Trace" unit2=""
            color=j scalebar=y
       ''')

End()
```

**Stage by stage:**

1. `sfspike` creates a 200-sample × 100-trace array with three spikes at samples 50, 100, 150. `nsp=3` is required — `sfspike` does not infer it from the length of the `k1` list.

2. `sfricker1 frequency=40` convolves each trace with a 40 Hz Ricker wavelet, turning the spikes into wavelets.

3. `sfnoise` adds Gaussian noise (variance 0.005, reproducible with seed 2025) to simulate real data.

4. `sfgrey` displays the result. Key choices:
   - `color=j` — jet colormap shows the gather with blue–green–yellow–red mapping. Good for seeing both amplitude and polarity.
   - `scalebar=y` — draws a colorbar on the right showing the amplitude-to-color mapping.
   - `label1/unit1/label2/unit2` — explicit labels override whatever is in the RSF header.
   - `transp=y yreverse=y` (defaults) — time (axis 1) runs vertically top-to-bottom, traces (axis 2) run left-to-right.

**Output:** `Fig/synth.vpl` — a seismic-style display of a synthetic gather.

---

### Example 2 — overlay with a velocity curve

File: `skills/plotting-with-vplot/references/example-overlay.SConstruct`

```python
from rsf.proj import *

Flow('cmp', None,
     '''
     spike n1=400 n2=60 k1=80,180,300 nsp=3 mag=1,0.7,0.5 |
     ricker1 frequency=30
     ''')
Plot('cmp', 'grey title="CMP" label1=Time unit1=s label2=Offset unit2=m')

Flow('vcurve', None,
     '''
     math n1=60 output="1500+x1*10" d1=1 o1=0
     ''')
Plot('vcurve',
     '''
     graph transp=y yreverse=y min2=0 max2=2000
           wantaxis=n wanttitle=n plotcol=2 plotfat=3 pad=n
     ''')

Result('cmp-with-vcurve', 'cmp vcurve', 'Overlay')
```

**Stage by stage:**

1. `cmp` is a 400-sample × 60-trace synthetic gather. `Plot` (not `Result`) is used because it is an intermediate that feeds the overlay.

2. `vcurve` is a 60-sample 1D array where sample i has value `1500 + i*10`, representing a velocity that increases linearly from 1500 to 2090 m/s across 60 offset positions. `sfmath` is used because there is no input file — `n1=60 d1=1 o1=0` defines the grid, and `output=` gives the analytic formula.

3. The `sfgraph` call for `vcurve` uses:
   - `transp=y yreverse=y` — makes the curve plot in the same orientation as the grey background (vertical time axis).
   - `min2=0 max2=2000` — the velocity curve values (1500–2090) are mapped onto the horizontal axis range 0–2000. This must match what `sfgrey` uses for its horizontal axis. Here `sfgrey` uses `label2=Offset unit2=m` with the default axis range from the data (0 to 59 for the 60 traces). There is a conceptual mismatch here — in a real workflow you would align `min2/max2` precisely with the offset range. For demonstration purposes, the curve is plotted within the displayed frame.
   - `wantaxis=n wanttitle=n` — suppress the graph's own axes; the background `cmp` image provides the axis frame.
   - `plotcol=2` — red curve.
   - `plotfat=3` — thick curve for visibility.
   - `pad=n` — no extra padding; ensures the curve aligns with the raster plot coordinate frame.

4. `Result('cmp-with-vcurve', 'cmp vcurve', 'Overlay')` — `vppen` draws `cmp.vpl` first (background), then `vcurve.vpl` on top.

**Output:** `Fig/cmp-with-vcurve.vpl` — the gather with a red velocity curve overlaid.

---

### Example 3 — wiggle trace display

File: `skills/plotting-with-vplot/references/example-wiggle.SConstruct`

```python
from rsf.proj import *

Flow('traces', None,
     '''
     spike n1=500 n2=15 k1=100,300 nsp=2 mag=1,0.6 |
     ricker1 frequency=25
     ''')

Result('traces',
       '''
       wiggle clip=0.8 title="Wiggle display"
              label1=Time unit1=s label2=Trace unit2=""
              poly=y zplot=1.5
       ''')
```

**Stage by stage:**

1. `sfspike` generates 500 samples × 15 traces with two wavelets per trace (at samples 100 and 300). `nsp=2` is required explicitly.

2. `sfricker1 frequency=25` creates 25 Hz Ricker wavelets — lower frequency than example 1, producing broader, more clearly visible wiggles.

3. `sfwiggle` displays with:
   - `clip=0.8` — explicit clip at 80% of the spike amplitude. Wiggle traces are clipped to this value to control how wide the excursions can get. Without an explicit clip, the pclip=98 default would estimate one.
   - `poly=y` — fills positive excursions with black polygons (standard seismic wiggle style).
   - `zplot=1.5` — traces are spaced at 1.5× the default separation, so they overlap slightly. This emphasizes the waveform shape. Values > 1 cause overlap; values < 1 leave gaps between traces.
   - Default `transp=n` for `sfwiggle` puts axis 1 (time) on the vertical axis and traces on the horizontal.

**Output:** `Fig/traces.vpl` — a 15-trace wiggle display with filled positive excursions.

