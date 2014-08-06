/** \mainpage IWAVE Trace I/O Package

Authors: William W. Symes and Tetyana Vdovina
<p>
This package uses interpolation and SU library functions (sucore) to
sample time traces from gridded spatial data and read/write these
traces to SEGY format files.  Main features:

- struct \ref tracegeom stores essential header and sampling info
for a single shot record, and a data buffer sized to accommodate
its data samples;

- function \ref init_tracegeom initializes a \ref tracegeom from a
file; may optionally initialize the data buffer from the file also.

- function \ref sampletraces extracts samples from a grid at
locations specified in a \ref tracegeom object, and stores the
samples in the object's data buffer.

- function \ref writetraces constructs SEGY traces from the data
in a \ref tracegeom and writes these to a specified file.

A "destructor" and "default constructor" are also provided, to give
\ref tracegeom and its associated functions as many of the
characteristics of a class as possible.

The central concept underlying the structure of this package is the
relation between GLOBAL and LOCAL grids. The GLOBAL grid describes a
(virtual) gridding of all of R^n; its characteristics are the
dimension, the grid steps along each axis, and the coordinates of the
grid origin (point with grid indices = 0). Sampling takes place from a
LOCAL grid, which is a (finite) subgrid of the GLOBAL grid. Data
corresponding to the LOCAL grid is assumed to reside in the memory
space of the process. The LOCAL grid must have the same grid steps
along each axis as the GLOBAL grid, and is additionally characterized
by the number of samples aong each axis and the coordinates of the
first sample. The function \ref init_tracegeom uses the descriptors of
GLOBAL and LOCAL grids together with trace headers read from an input
file to compute indices and relative in-cell coordinates of each
receiver point and of the source point. It is assumed that the input
file describes a single shot record, with a unique source location,
and this is read from the first trace in the input header file. This
sampling information is stored in a \ref tracegeom object, which is
subsequently used by other functions in the package to extract samples
from a (LOCAL) grid of data, and to construct trace headers for trace
output.

The time sample rate for simulation is an independently prescribed
input - it is presumed that a sample rate is set somehow in the setup
phase of the simulator, and is available to the trace package as an
input. Sampling to or from the internal trace buffer of \ref tracegeom
occurs at this simulation sample rate. The function \ref writetraces
uses cubic spline interpolation (see \ref cubic) to resample this data
as prescribed by the trace headers used to initialize the \ref
tracegeom object. Thus output traces with exactly the same geometry as
these input headers. An important exception: only those traces with
receiver positions lying inside of the LOCAL grid will be sampled -
traces for receiver positions outside of the LOCAL grid are neither
sampled, nor read from an input file, nor written to an output file.

Alternatively, the user may specify the number of time steps and first
sample time, in which case the external traces have this length and are
sampled at the internal simulation rate. This option is use when
precise control of simulation time step is essential, for instance in
some convergence studies. See the documentation of \ref init_tracegeom
for precise details.

Parallel implementation: all i/o takes place on rank 0. Reads: data
read trace-by-trace into segy trace buffer from file on rank 0, then
broadcast to all other ranks. Each rank computes that part of the data
needed locally, then stores it. Writes: roughly the reverse of read,
but implemented with point-to-point communication, hence slower.

Source and receiver coordinates may be interpreted as relative to a 
coordinate offset vector, an argument to the construct function. This
feature is useful in using this data structure to implement a towed
streamer source, in which the input coordinate vector is interpreted 
as the coordinate of a reference point on the source array, and the
receiver coordinates of the data as relative coordinates to it. Thus 
for example a superposition of three point sources, at 10 m to left and 
right of a central source, might be implemented by (1) creating a SEGY
data file with three traces, having gx=-10, 0, and 10 respectively, and
(2) passing the coordinates of the source reference point via the auxiary 
offset vector. Alternatively, the same effect can be achieved by passing 
a zero offset vector (for example an RPNT_0) and adding the reference
point coordinates to the trace header gx.

Note above all that, in defining sources using the load option with the 
traceio package, the source coordinates play no role at all - the 
receiver coordinates give the source positions, at which data is added 
into the grid, either relative to a source reference point or absolute 
in global grid coordinates. 

<hr>
<a href="../../../doc/html/index.html">IWAVE Home Page</a>
*/

