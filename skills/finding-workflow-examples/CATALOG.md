# Workflow & Dataset Catalog

Curated index into `book/` — the Madagascar corpus of runnable SConstructs. This catalog names one or two canonical examples per domain; it is not exhaustive. For anything not listed, grep `book/` directly (see `SKILL.md`).

Paths are relative to the repo root (`/Users/jgoai/m8r/src/`). Open the file itself for the authoritative flow; the one-line descriptions here are pointers only.

## Workflows

### Synthetic modeling

- **Tuning wedge** — `book/rsf/tutorials/wedge/SConstruct` — variable-thickness wedge model for resolution analysis.
- **Ricker wavelet analysis** — `book/rsf/tutorials/wavelet/SConstruct` — build and characterize a Ricker wavelet, frequency content.
- **Finite-difference acoustic modeling (BP 2004)** — `book/rsf/school/modeling/SConstruct` — 2D acoustic wave propagation on the BP 2004 velocity model using awefd2d.
- **Simple layered / impulse model** — `book/rsf/school/marm/SConstruct` — eikonal traveltime and anisotropic moveout on the Marmousi velocity model.

### NMO & velocity analysis

- **NMO starter** — `book/rsf/tutorials/nmo/SConstruct` — apply NMO to a synthetic CMP using sfnmo, step-by-step.
- **Velocity scan + semblance + stack** — `book/rsf/tutorial2017/field/SConstruct` — vscan, semblance picking, NMO correction, and hyperbolic Radon demultiple on a Viking Graben CMP.

### Migration

- **Marmousi migration (one-way wave-equation)** — `book/gallery/marmousi/oway/SConstruct` — extended split-step (zomig3) migration on the Marmousi exploding-reflector dataset.
- **Sigsbee migration (school)** — `book/rsf/school/sigsbee/SConstruct` — acoustic finite-difference modeling on the Sigsbee 2A model.
- **Reverse-time migration helper** — `book/Recipes/rtm.py` — RTM building blocks called from SConstructs.
- **One-way wave-equation migration helper** — `book/Recipes/wemig.py` — oneway WE migration.
- **Zero-offset migration helper** — `book/Recipes/zomig.py` — zero-offset migration.
- **Pre-stack migration helper** — `book/Recipes/pmig.py` — prestack migration.

### Denoising & filtering

- **Bandpass filtering (inline)** — `book/rsf/school2020/seismic/SConstruct` — sfbandpass applied as part of a 2D land-data processing flow.
- **F-x deconvolution** — `book/zju/optnoise/postack/SConstruct` — post-stack noise attenuation using fxdecon together with dip-guided EMD.
- **Dip filtering** — `book/rsf/su/rsfdipfilt/SConstruct` — noise attenuation by F-K dip filter.

### Transforms

- **Hyperbolic Radon transform** — `book/rsf/tutorial2017/radon/SConstruct` — forward/adjoint hyperbolic Radon.
- **2D FFT / F-K (inline)** — `book/rsf/tutorials/timefreq/SConstruct` — fft1, spectra, and F-K filtering demonstrated on sinusoid test signals.
- **Hilbert transform / envelope / instantaneous attributes** — `book/rsf/tutorials/hilbert/SConstruct`.
- **Time-frequency analysis** — `book/rsf/tutorials/timefreq/SConstruct` — S-transform (st), STFT, and local time-frequency decomposition (timefreq program).

### Interpolation

- **PEF-based interpolation** — `book/rsf/tutorials/interp/SConstruct`.
- **Spitz interpolation** — `book/rsf/tutorials/spitz/SConstruct`.

### Well-tie / log calibration

- **Well-to-seismic tie** — `book/rsf/tutorials/well-tie/SConstruct`.

### Geometry & survey setup

- **Regular 2D survey** — `book/rsf/tutorials/survey/SConstruct` — define source/receiver grid.
- **Real SPS/XPS/RPS geometry** — `book/rsf/school2020/seismic/SConstruct` — parse real-world survey geometry files into RSF.

### Rock physics & petrophysics

- **Rock physics relations** — `book/rsf/tutorials/rockphysic/SConstruct`.
- **Log conditioning (v1)** — `book/rsf/tutorials/petro1/SConstruct`.
- **Log conditioning (v2)** — `book/rsf/tutorials/petro2/SConstruct`.

### Facies, attributes, ML

- **Facies classification** — `book/rsf/tutorials/facies/SConstruct`.
- **Seismic attributes** — `book/rsf/tutorials/attr/SConstruct`.
- **Neural network on seismic** — `book/rsf/tutorials/nn/SConstruct`.

### Plotting idioms

- **Colormaps (including scientific palettes)** — `book/rsf/tutorials/colormaps/SConstruct`.
- **Colored sections** — `book/rsf/tutorials/colored/SConstruct`.
- **Lineament tracking on images** — `book/rsf/tutorials/lineaments/SConstruct`.
