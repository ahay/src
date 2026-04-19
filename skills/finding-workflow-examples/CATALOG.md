# Workflow & Dataset Catalog

Curated index into `book/` — the Madagascar corpus of runnable SConstructs. This catalog names one or two canonical examples per domain; it is not exhaustive. For anything not listed, grep `book/` directly (see `SKILL.md`).

Paths are relative to the repo root (`/Users/jgoai/m8r/src/`). Open the file itself for the authoritative flow; the one-line descriptions here are pointers only.

## Workflows

### Synthetic modeling

- **Tuning wedge** — `book/rsf/tutorials/wedge/SConstruct` — variable-thickness wedge model for resolution analysis.
- **Ricker wavelet analysis** — `book/rsf/tutorials/wavelet/SConstruct` — build and characterize a Ricker wavelet, frequency content.
- **Finite-difference acoustic modeling (BP 2004)** — `book/rsf/school/modeling/SConstruct` — 2D acoustic wave propagation on the BP 2004 velocity model using awefd2d.
- **Analytic velocity model + ray tracing** — `book/rsf/school/ray/SConstruct` — gradient velocity model with a Gaussian anomaly; demonstrates rays2 ray-tracing on a 2D field.

### NMO & velocity analysis

- **NMO starter** — `book/rsf/tutorials/nmo/SConstruct` — apply NMO to a synthetic CMP using sfnmo, step-by-step.
- **Velocity scan + semblance + NMO** — `book/rsf/tutorial2017/field/SConstruct` — vscan, semblance picking, NMO correction, and hyperbolic Radon demultiple on a Viking Graben CMP.

### Migration

- **Marmousi migration (one-way wave-equation)** — `book/gallery/marmousi/oway/SConstruct` — extended split-step (zomig3) migration on the Marmousi exploding-reflector dataset.
- **Sigsbee migration (one-way wave-equation)** — `book/gallery/sigsbee/oway/SConstruct` — split-step (zomig3) one-way wave-equation migration on the Sigsbee zero-offset dataset.

The following four entries are importable Python helpers (not standalone SConstructs) — see SKILL.md adaptation step 4 for how to use them.

- **Reverse-time migration helper** — `book/Recipes/rtm.py` — RTM building blocks called from SConstructs.
- **One-way wave-equation migration helper** — `book/Recipes/wemig.py` — oneway WE migration.
- **Zero-offset migration helper** — `book/Recipes/zomig.py` — zero-offset migration.
- **Pre-stack migration helper** — `book/Recipes/pmig.py` — prestack migration.

### Denoising & filtering

- **F-x deconvolution** — `book/zju/optnoise/postack/SConstruct` — post-stack noise attenuation using fxdecon together with dip-guided EMD.
- **Dip filtering** — `book/rsf/su/rsfdipfilt/SConstruct` — noise attenuation by F-K dip filter.

### Transforms

- **Hyperbolic Radon transform** — `book/rsf/tutorial2017/radon/SConstruct` — forward/adjoint hyperbolic Radon.
- **Hilbert transform / envelope / instantaneous attributes** — `book/rsf/tutorials/hilbert/SConstruct` — uses envelope with hilb=y to compute Hilbert transform, envelope, and instantaneous phase on a seismic trace.
- **Time-frequency analysis** — `book/rsf/tutorials/timefreq/SConstruct` — S-transform (st), STFT, and local time-frequency decomposition (timefreq program).

### Interpolation

- **PEF-based interpolation** — `book/rsf/tutorials/interp/SConstruct` — compares linear (remap1), sinc-lag (inttest1), and spline 1D interpolation methods.
- **Spitz interpolation** — `book/rsf/tutorials/spitz/SConstruct` — uses hpef and signoi to separate signal and noise via helix PEF estimation.

### Well-tie / log calibration

- **Well-to-seismic tie** — `book/rsf/tutorials/well-tie/SConstruct` — reads LAS file via las2rsf, despiked DT/RHOB logs, and computes a synthetic seismogram.

### Geometry & survey setup

- **Regular 2D survey** — `book/rsf/tutorials/survey/SConstruct` — define source/receiver grid.
- **Real SPS/XPS/RPS geometry** — `book/rsf/school2020/seismic/SConstruct` — parse real-world survey geometry files into RSF.

### Rock physics & petrophysics

- **Rock physics relations** — `book/rsf/tutorials/rockphysic/SConstruct` — classifies sand/shale facies from well CSV and plots IP/VP/VS cross-plots.
- **Log conditioning (v1)** — `book/rsf/tutorials/petro1/SConstruct` — classifies brine-sand/oil-sand/shale from VSH/SW/PHI/IP well logs and plots histograms.
- **Log conditioning (v2)** — `book/rsf/tutorials/petro2/SConstruct` — extends petro1 with Gassmann fluid substitution (brine-to-oil) on VP/VS/RHO logs.

### Facies, attributes, ML

- **Facies classification** — `book/rsf/tutorials/facies/SConstruct` — loads wireline facies-label CSV (csv2rsf); minimal entry point for facies ML workflows.
- **Seismic attributes** — `book/rsf/tutorials/attr/SConstruct` — computes RMS, absolute-mean, and envelope-mean amplitude attributes using envelope and sfattr.
- **Neural network on seismic** — `book/rsf/tutorials/nn/SConstruct` — plots sigmoid/tanh/ReLU activation functions and applies a neural network to well-log data (las2rsf).

### Plotting idioms

- **Colormaps (including scientific palettes)** — `book/rsf/tutorials/colormaps/SConstruct` — demonstrates spectral, linearlfb, and viridis palettes using grey on random data.
- **Colored sections** — `book/rsf/tutorials/colored/SConstruct` — fetches an acoustic-impedance section and well log; displays with inferno colormap via grey.
- **Lineament tracking on images** — `book/rsf/tutorials/lineaments/SConstruct` — computes directional derivatives of a Bouguer gravity map to extract structural lineaments.
