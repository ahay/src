/** \mainpage IWAVE Isotropic Elastic Wave Equation Solver
<p> 
Purpose: solve the isotropic wave equation in 2 or
3 spatial dimensions, output stress traces at specified sample rates
and/or movie frames of pressure. Uses staggered grid finite difference
scheme of order 2 in time and 2k in space, k=2, derived from
stress-velocity form of elastic dynamics. Either pressure-free
(reflecting) or absorbing (PML) boundary conditions on boundary faces
of simulation hypercube.
<p>
Authors: Igor S. Terentyev, Tetytana Vdovina, William W. Symes, Xin Wang
<p>
Usage: 
  <ul>
  <li>esg.x datafile= [optional parameters], or </li>
  <li>esg.x par=[name of parameter file]</li>
  </ul>

Due to the number of parameters, the second form is usually the more
convenient, that is, store all of the parameters in a file. The
command name (esg.x) is a standin for any valid form of invocation,
such as MPI-enabled invocation using mpirun or mpiexec, with
appropriate path information. 

<p>

Table of Contents:
<p>
\ref typical
<p>
\ref alternate
<p>
\ref esgnotes
<p>
<p>

<hr>
<a href="../../../doc/html/index.html">IWAVE Home Page</a>

*/
