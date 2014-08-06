/** \mainpage IWAVE Trace Sampler Package

Authors: William W. Symes and Tetyana Vdovina
<p>
This package defines the interactions between gridded spatial data
and trace time series.

The TRACE_TERM struct stores a \ref tracegeom (trace header data 
structure) set up to set up to sample a specific field (local gridded
data). Both local and global grids are attributes of an IMODEL; the 
necessary data is extracted in the \ref traceterm_init function and used 
to construct the captive \ref tracegeom. Data is either inserted in grid 
(load=1) or extracted from grid (load=0), via call to \ref sampletraces.
The \ref traceterm_init function also permits the user to offset all 
source and receiver coordinates by a fixed vector, which is useful in
defining towed streamer sources.

SAMPLER is a shallow wrapper around TRACE_TERM. Besides an internal
TRACE_TERM, SAMPLER also stores an array (\ref IPNT) of field indices,
The run method cycles through the fields indexed by the elements of
this array. This feature is essential in defining source behaviour (load 
option) for systems such as global PML which duplicate dynamical fields 
in the physical domain. For such systems, source updates must be made to 
all copies of the relevant fields. For receiver behaviour (non-load or save
option), only the first (0) indexed field is sampled.

MOVIE defines a sequence of snapshots of dynamic fields. Choice of
fields to sample, and mnemonic keys used to specify this choice,
deferred to specific model packages, but a simple interface for this
key definition is provided. Thus MODEL plays the role of an abstract
base.

<p>
Table of Contents:
<p>
\ref trace_term
<p>
\ref sampler
<p>
\ref movie

<hr>
<a href="../../../doc/html/index.html">IWAVE Home Page</a>

*/

