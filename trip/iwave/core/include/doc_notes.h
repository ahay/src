/**
@page notes Notes on Modeling with IWAVE 
<p> 
<ul> 

<li>Survey Simulation: IWAVE modeling applications simulate all shots
defined by the source-receiver geometry of a survey, either serially
or in parallel. To simulate multiple shots, simply supply an output
file, identified by the datafile parameter, or a prototype header
file, containing SEGY trace headers specifying the survey geometry,
sorted in shot order.  Each shot will be simulated and overwritten on
the output file traces in the order of appearance in the file. In
serial mode, shots are simulated one at a time. In distributed (MPI)
mode, as many shots are simulated simultaneously as possible, for the
number of processes assigned to the job. See \ref parallel for more on
this. </li>

<li> Parameter Parsing: Parameters are specified as key=value pairs, as in SU and
similar applications. String parameter values representing file names
can be either names of files in the working directory, or relative or
fully qualified pathnames. Any line not containing the separator
character (by default, the equals sign =) is treated as a
comment. Whitespace is ignored, including whitespace between key and
value strings and the separator (=). Key and value strings follow
regular expression syntax - in particular, key or value strings with
embedded whitespace must be quoted. Inline text separated from key or
value strings by whitespace is treated as comment. Lines can be
commented out by quoting them. The parameter files generated in the 
<a href="../../../demo/doc/html/index.html">IWAVE demo package</a> 
illustrate some of these features. For further information about
parameter parsing in IWAVE, see 
<a href="../../../base/doc/html/index.html">IWAVE base package documentation</a> .
</li> 
<p>
<li> 
Internal units - chosen to conform to unit
choices made, explicitly or implicitly, in the SEGY standard. <ul>
<li>length = m</li> <li>time = ms</li> <li>mass = kg</li> </ul> Note
that frequency should be specified in KHz = cycles/ms, rather than
Hz. Other consistent derived units and typical ranges for sedimentary
rocks are <ul> <li>velocities = O(1e+0) m/ms,</li> <li>bulk moduli =
O(1.e+3 - 1.e+4) MPa,</li> <li>densities = O(1.e+3) kg/m^3.</li> </ul>
All data must be scaled to these units before or during input. Note
that an IWAVE extension to RSF allows the user to change metric units
'on the fly', thus avoiding the need to preprocess data files. To
scale RSF data during input, use the scale keyword to specify a power
of 10 scaling factor. For example, if velocity samples are in m/s, add
'scale=-3' to the RSF header file for velocity, so that the data is
recorded internally in the expected m/ms units.  </li> <p> 

<li>
Axes: Numbers (rather than letters) describe axes. We use the RSF/SEP 
system numbering (starting with 1). The identity of the axes (which one 
is depth, etc.) is vital to establish the correct sampling of spatial 
fields to SEGY trace output. IWAVE's default axis assignment is 1=z, 
(fastest axis), 2=x, 3=y. However, another IWAVE extension to RSF allows 
the user to override these defaults. If for example the axis ordering for 
density is 1=x (fastest axis), 2=y, and 3=z (as is conventional for GOCAD 
and many other applications), simply add "x_axis=1", "y_axis=2", and 
"z_axis=3" to the RSF header file. This facility allows the user to avoid
transposing data files.
<p>
Also note that only the information corresponding to active axes need be 
present in RSF header files. 1D simulation requires only the specification 
of axis 1 quantities,
2D only axes 1 and 2.
  </li>
<p>
<li>
Word order: IWAVE uses the RSF data_format keyword to control input word order
for material parameter data. Two choices are provided: native_float and xdr_float.
Since XDR word order is IEEE big endian, the choice data_format=xdr_float permits
big endian data to be read on little endian machines. The cost of the XDR stream
operations, and the proportion of cycles spent reading material parameter data,
are both so small in a typical 3D simulation that we strongly recommend that all 
material parameter data be XDR encoded.
<p>
Likewise, <a href="../../../trace/doc/html/index.html"> IWAVE's trace i/o package</a> 
uses SU i/o functions, hence inherits SU's choice of XDR or native floating point word
order. This choice is made at compile time, rather than run time as with RSF, via the SUXDR
compiler directive (see the <a href="../../../README.INSTALL">IWAVE installation notes</a> 
for more information).
  </li>

<p>
 
  <li>
In all IWAVE apps, verbose diagnostic output (switched by the dump flags) is directed to a file 
named "cout[rk].txt" in the working directory, where [rk] is a fixed-width 
representation of process rank. Thus in serial execution verbose output 
   goes to cout0.txt. In parallel (domain-decomp) execution with 2048 processes, root output goes to cout0000.txt, process 123 output goes to cout0123.txt, etc. The software creates the filenames and files.

    <ul>
      <li>
      dump_pi   = [int] [Default = 0] - dump parallel/dom decomp info</li>
      <li>
      dump_lda  = [int] [Default = 0] - dump grid data for allocated arrays, 
      including ghost cell zones.</li>
      <li>
      dump_ldc  = [int] [Default = 0] - dump grid data for arrays updated 
      in timesteping (computed subarrays, not including ghost cells).</li>
      <li>
      dump_lds  = [int] [Default = 0] - dump grid data for send arrays</li>
      <li>
      dump_ldr  = [int] [Default = 0] - dump grid data for receive arrays</li>
      <li>
      dump_term = [int] [Default = 0] - dump sampling data (source, trace,
      movie, etc.)</li>
      <li>
      printact  = [int] [Default = 0] - timestep action output: 
        <ul>
	<li> 0 - no output </li>
        <li> 1 - timestep number</li>
        <li> 2 - detailed actions within timestep.</li>
        <li> 3 - all allocated arrays printed at every action - do this
        only for VERY small examples!!!</li>
        </ul>
      </li>
    </ul>
  </li>
  </ul>
*/
