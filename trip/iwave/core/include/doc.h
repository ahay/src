/** \mainpage IWAVE Simulator data structure package

Author: Igor Terentyev<br>
Other Contributors: Tetyana Vdovina, William W. Symes, Xin Wang
<p>
This package defines the data structure \ref IWAVE for a timestepping simulator on a regular spatial grid.
\ref IWAVE consists pricipally of two structs: \ref IMODEL, a coherent collection of \ref RDOMAIN structs defining arrays and virtual subarrays participating in the simulation, and \ref FD_MODEL, containing function pointers referenced in the functions implementing the simulation.

To create a simulator using IWAVE, it is necessary to define functions conforming to the interfaces listed in \ref FD_MODEL and implementing the computations described in its documentation. These include time step and i/o functions and various auxiliary functions which provide the basic information necessary to set up ghost cell arrays for boundary conditions and data exchange in domain decomposition.
<p>
\ref parallel gives operational details for parallel simulation via domain decomposition and simultaneous shots.
<p>
\ref notes overview parameter parsing, units, word order, and other details.
<p>
Several other packages support this one, including
<ul>
<li> 
<a href="../../../base/doc/html/index.html">base</a> - parameter parsinng, file handling, and other utilities
</li>
<li>
<a href="../../../grid/doc/html/index.html">grid</a> - grid (RSF) data structures and i/o functions 
</li>
<li>
<a href="../../../trace/doc/html/index.html">trace</a> - trace (SEGY) data structures and i/o functions
</li>
<li>
<a href="../../../sample/doc/html/index.html">sample</a> - functions for sampling trace data from/to grid data
</li>
</ul>
<hr>
<a href="../../../doc/html/index.html">IWAVE Home Page</a>

*/
