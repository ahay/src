/**
@page alternate Alternate Use Cases

<p>
This section describes various alternative tasks implemented by esg.x, and 
corresponding parameter selection patterns.
<ol>
<li>Create/modify output data file: it is possible to take headers from one 
file of SEGY traces, yet write the output to another. To do this, specify 
the source of the trace headers as hdrfile, and the output data file as 
datafile (as in the typical use case):
<ul>
<li>hdrfile = [string]</li>
<li>datafile = [string]</li>
</ul>
In this case, the output file is overwritten in its entirety (including 
any SEGY header information).
</li>
<p>
<li>
Movie output - uses IWAVE movie output facility, see <a
href="../../../sample/doc/html/index.html">IWAVE sampler package
documentation</a> for more on the structure of this facility. 

Usage: 
<ol>
<li>add "movie[n]=[code]" lines to par file to specify field(s) to be sampled. Legit codes for ASG model: p (pressure), v1, v2, v3 (components of velocity, in the internal axis ordering, i.e. v1=z-component, v2=x-component, v3=y-component. Any other code will be ignored. Thus for example
<p>
movie1=p
movie2=v3
<p>
will produce two movies, one of pressure, one of the y-component of velocity.</li>
<li>add "moviestep=[real]" line to par file to specify time step between movie frames, in units of ms. Thus to sample ten frames per second, include "moviestep=100.0". Default = DEFSTEP*dt, where dt is the <i>internal</i> time step of the simulartor, and DEFSTEP is defined in <a href="../../../sample/include/movie.h">sample/include/movie.h</a> (so this default is not very useful, since it depends on dt which one doesn't know a priori - better to specify moviestep).</li>
<li>only 2D movies are produced. Therefore no 1D movies at all, and only a single plane may be sampled in the 3D case. The sampling plane is perpindicular to one of the coordinate axes; which axis it's perp to is the value of the "movieaxis3d=[int]", legit values 0,1,2, default = 2 (this choice refers to IWAVE's internal axis ordering, so in the standard (z,x,y) ordering it's the (z,x) plane, but in the other common ordering (x,y,z) it's the (x,y) plane - depends on the x-axis, y-axis, and z-axis keyword settings, see <a href="../../../grid/doc/html/index.html">IWAVE Grid i/o package docs</a>).</li>
<li>Also specify the position of the sampling plane by where it cuts the perpindicular coordinate axis, by "movieslice3d=[real]". Default = plane with minimum perp coordinate, i.e. a face of the simulation cube.</li>
</ol>
<p>
Each movie is written to an RSF file named movie_[code].rsf(@). The RSF header 
contains the information required to display the movie using any appropriate utility.
For instance, a movie of pressure snapshots with header file movie_p.rsf might contain
n1=101 n2=201, in which case one could view the corresponding movie with
<p>
xmovie n1=101 n2=201 [other params] < movie_p.rsf@
</li>
<p>
<li>
Alternate methods for specifying computational grid and variable compressional velocity and shear velocity. 
It is possible to prescribe it directly:
<ul>
<li>pvelocity = [string], and svelocity = [string]</li>
</ul>
</li>
<p>
<li>Alternate method for specifying computational grid and constant material 
parameters:  specifies grid in fashion of RSF header, constant values for 
bulk modulus and bouyancy (or equivalent parameters):
<ul>
<li>n1 = [int] [default = 1],</li>
<li>n2 = [int] [default = 1],</li>
<li>n3 = [int] [default = 1],</li>
<li>d1 = [float] [default = 1.0],</li>
<li>d2 = [float] [default = 1.0],</li>
<li>d3 = [float] [default = 1.0],</li>
<li>o1 = [float] [default = 0.0],</li>
<li>o2 = [float] [default = 0.0],</li>
<li>o3 = [float] [default = 0.0],</li>
<li>x_axis = [int] [default = 2],</li>
<li>y_axis = [int] [default = 3],</li>
<li>z_axis = [int] [default = 1],</li>
<li>reflambda = [float] [default = 2250 MPa],</li>
<li>refmu    = [float] [default = 2250 MPa],</li>
<li>refbouy  = [float] [default = 0.001 m^3/kg],</li>
</ul>
Notes: axis information supplied mostly for test purposes - since the material 
modeled is homogeneous, the axis labeling is immaterial. The values of Lame's parameters 
(reflambda and refmu) and bouyancy (refbouy) have priority: the other two are used
only if one of these is not given.
</li>
<p>
<li>Direct specification of bouyancy from RSF file, rather than density:
<ul>
<li>bouyancy = [string]</li>
</ul>
</li>
<p>
<li>Direct specification of bouyancy or density on shifted grids:
<ul>
  <li>bouyancy1 = [string] and</li>
  <li>bouyancy2 = [string] and</li>
  <li>bouyancy3 = [string], </li>
  </ul>
  or
  <ul>
  <li>density1  = [string] and</li>
  <li>density2  = [string] and</li>
  <li>density3  = [string] </li>
  </ul>
  or
  <ul>
  <li>rho1      = [string] and</li>
  <li>rho2      = [string] and</li>
  <li>rho3      = [string] </li>
  </ul>
For this option, the RSF files should represent shifted values of the fields, by 
one-half grid cell along each axis. The option is provided mostly for test purposes:
it allows the use of directly sampled, rather than interpolated, values of bouyancy
in case the latter is somehow defined as an actual function. Thus interpolation 
does not intervene in the convergence behaviour of the scheme.
</li>
<p>

<li>Right hand side of wave equation read directly from file: the
typical use case specifies not the time function (source pulse)
appearing in the point defect model of an isotropic point radiator
(multiplying a spatial delta on the RHS of the the constitutive law,
relating rate of change of pressure to rate of change of volume), but
rather the pulse produced at a given distance from the source point in
a homogeneous medium. This is a convenient and meaningful calibration,
but for some purposes it may be more appropriate to specify the source
pulse directly. To enable this option, simply specify
<ul>
<li>refdist = 0.0</li>
</ul>
with other inputs as in the typical use case.
</li>
<p>
<li>Ricker pulse with given peak frequency: to produce a Ricker pulse
(in 3D), of specified peak frequency and unit amplitude at specified distance
from source point: inputs as in source specification in typical use
   case, with these amendments: 
<ul>
<li>omit source keyword</li>
<li>fpeak = [float] [default = 0.01 KHz] (also included in PML inputs) 
</ul>
The reference distance specified by keyword refdist must be positive in this
case. The pulse produced in 2D is the 3/2-derivative of a Gaussian, rather
than a Ricker. 
</li>
<p>
<li>Array source: source is an array of point sources, implicitly summed to produce RHS input to the pressure wave equation. Sources act at the prescribed points on the RHS, with adjoint interpolation to gridpoints if sampord=1.
<ul>
<li>srctype=array</li>
<li>source=[filename]</li>
</ul>
The file holding the array source data should be an SU (SEGY trace) file. The traces should have their (source) locations identified by <i>receiver</i> coordinates (that is, keywords gx, gy, gelev) rather than source coordinates (sx, sy, selev). For array source files, the source keywords have no significance. The reason for this choice is that arrays sources are used in RTM applications of IWAVE, in which case receiver locations act as sources. Choosing the receiver header words to hold array source information avoids an unnecessary data transformation step in the RTM application, and is just as convenient as the other obvious option in general.
<p>
Array sources may be used to simulation incident plane waves or other extended source effects, or multipole near-point sources with nontrivial radiation patterns.
<p>
An array trace with one trace produces exactly the same data as a point source with the same trace input and refdist=0.0.
</li>
<p>
<li>
Alternate methods for specifying time steps: see \ref esgnotes.
</li>
</ol>
*/
