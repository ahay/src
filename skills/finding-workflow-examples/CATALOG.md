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

## Datasets

Two styles of dataset recipe live in `book/`:

- **Gallery Python helpers** (`book/gallery/<name>/<name>.py`) expose functions like `getvel('vel_target')`, `get_shots('shots_target')`, `get_zodata(...)` that internally call `Fetch()` and `Flow()`. Use them by adding the gallery directory to `sys.path` and `import <name>` in your SConstruct.
- **`book/data/<name>/` directories** contain fetch-and-prep SConstructs for real or published datasets. Read the SConstruct for the `Fetch()` server and the convert-to-RSF chain.

### Gallery helpers

- **Marmousi** — `book/gallery/marmousi/marmousi.py` — `getvel`, `get_zodata`, `get_shots`, `get_ffd_shots`; velocity, zero-offset exploding-reflector data, and prestack shot records.
- **Sigsbee 2A** — `book/gallery/sigsbee/sigsbee.py` — `getvel(vel, veltype)` (migvel or strvel), `getzo`, `gethrzo`, `getshots`; canonical deep-water salt-body benchmark.
- **BP 2004 velocity benchmark** — `book/gallery/bp/bp.py` — `getvel`, `getden`, `getshots`, `zodata`; velocity and density grids plus prestack shot records for the BP 2004 tomography benchmark.
- **BP TTI 2007** — `book/gallery/bptti/bptti.py` — `getmod(par)` (par in epsilon/delta/vp/theta), `getshots`; anisotropic (TTI) parameter grids and prestack shots for the BP 2007 anisotropy benchmark.
- **Pluto 1.5** — `book/gallery/pluto/model/SConstruct` — no Python helper; fetch-and-convert SConstruct for the Pluto 1.5 velocity model (depth-interval velocity SEGY from the pluto server).
- **SEG/EAGE Overthrust** — `book/gallery/overthrust/overthrust.py` — `getvel`; 3D thrust-belt P-wave velocity cube fetched from the SEG/EAGE open-data S3 bucket.
- **Teapot Dome 3D** — `book/gallery/teapot/teapot.py` — `get_vrms1`, `get_vint1`; RMS and interval velocity functions for the Teapot Dome 3D survey geometry benchmark.
- **French dome model** — `book/gallery/french/french.py` — `get_refl`, `get_zodata2d`, `get_zodata`; synthetic French dome reflectivity surface and corresponding zero-offset and prestack data.
- **Hess VTI** — `book/gallery/hessvti/hessvti.py` — `get_model`, `get_shots`, `get_zodata`; builds vp, delta, epsilon, crho grids at import time; VTI anisotropy benchmark.
- **SEG/EAGE Salt** — `book/gallery/segsalt/segsalt.py` — `getvel2D`, `getvel3D`; 2D and 3D velocity fields for the SEG/EAGE salt model fetched from the open-data S3 bucket.
- **Two-half-space model** — `book/gallery/twohalf/twohalf.py` — `getvel`, `getshots`; simple two-half-space synthetic velocity model and prestack shot records.
- **v(z) layered model** — `book/gallery/vofz/vofz.py` — `get_velocity`, `zero_offset`, `shots`, `cmps`, `get_impulse`, `impulse_response`; v(z) gradient layered velocity with layered synthetic data generation and migration helpers.
- **Constant-velocity model** — `book/gallery/constant/constant.py` — `get_zodata`, `get_cmps`, `get_impulse`; constant-velocity synthetic dataset used for impulse-response and migration tests.
- **Gradient-velocity model** — `book/gallery/gradient/gradient.py` — `get_velocity`, `zero_offset`, `shots`, `cmps`, `get_impulse`, `impulse_response`; linear-gradient velocity model with synthetic data and migration helpers, analogous to vofz.
- **1994 Amoco model** — `book/gallery/model94/model94.py` — `get_vel`, `get_shots`, `get_shot_headers`; the 1994 BP/Amoco velocity model and prestack shot records from the bpmodel94 open-data archive.
- **1994 statics challenge** — `book/gallery/statics94/SConstruct` — no Python helper; gallery-level SConstruct for the 1994 BP statics benchmark (see also `book/data/bpstatics94/` for the fetch-and-prep recipes).

### data/ fetch-and-prep recipes

- **Marmousi 2 (elastic)** — `book/data/marmousi2/` — elastic version of Marmousi with density and Vs; includes `fdMod`, `vx`, `vz`, and `div` subdirectories for elastic finite-difference modeling.
- **SEAM Phase 1 2D** — `book/data/seam-phase1-2d/` — synthetic SEAM Phase 1 2D classic slice; includes fetch, filter, forward modeling, and vscan subdirectories.
- **New Zealand Kahu 3D** — `book/data/new-zealand-kahu-3d/` — representative of New Zealand 3D surveys; siblings `new-zealand-kerry3d`, `new-zealand-opunake-3d`, `new-zealand-parihaka-3d`, `new-zealand-tui3d`, `new-zealand-waihapa-3d`, `new-zealand-waipuku-3d`, `new-zealand-waka-3d` follow the same pattern.
- **Chevron 2013 benchmark** — `book/data/chevron2013/` — Chevron 2013 FWI benchmark dataset with fetch and firstlook subdirectories.
- **Chevron 2014 benchmark** — `book/data/chevron2014/` — Chevron 2014 benchmark with fetch, fwi, geom, and survey-volume subdirectories.
- **FreeUSP 2D Land data** — `book/data/freeusp/` — used by `book/rsf/school2020/seismic/SConstruct`; fetch and land subdirectories for the FreeUSP land seismic dataset.
- **Alaska OBS data** — `book/data/alaska/` — ocean-bottom seismometer data from Alaska; contains SConstruct and gathers subdirectory.
- **Oz land data** — `book/data/oz/` — Australian land seismic lines (line31-81, line41-81 variants) with SU-to-RSF conversion and NMO processing examples.
