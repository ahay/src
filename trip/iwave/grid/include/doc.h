/** \mainpage IWAVE Grid I/O package

Authors: William W. Symes and Xin Wang
<p>
Implements regular grid data structure and related i/o, based on RSF 
file format (http://www.reproducibility.org/wiki/Guide_to_RSF_file_format)
<p>
<ul>
<li> struct \ref axis defines a \ref grid axis. An \ref axis has four attributes: the
number of gridpoints on the axis (n), the step between gridpoints (d),
the coordinate of the first gridpoint (o), and the axis index (id)
which indicates its position in the a reference ordering of
axes (see Notes below for more on this). Functions provided:
<ul> 
<li>\ref init_default_axis : default initialization </li>
<li>\ref fprint_axis< print axis to stream /li>
<li>\ref print_axis : print axis to stdout </li>
<li>\ref compare_axis : compare two axes for equality </li>
</ul>
<li> struct \ref grid consists of up to \ref IWAVEBASE::RARR_MAX_NDIM axes. Functions 
provided:
<ul>
<li>\ref init_grid : initializes grid of dimension dim - default-initializes dim axes, assigns id's in natural order </li>
<li>\ref fprint_grid : print grid to stream </li>
<li>\ref fprint_grid : print grid to stdout </li>
<li>\ref compare_grid : compare two grids by comparing dimensions, axes
<li>\ref get_datasize_grid : return number of gridpoints </li>
<li>\ref get_n : extract array (\ref IPNT) of axis lengths (n's) </li>
<li>\ref get_d : extract array (\ref RPNT) of grid steps (d's) </li>
<li>\ref get_o : extract array (\ref RPNT) of grid origin coordinates (o's) </li>
<li>\ref get_gs : extract array (\ref IPNT) of global indices of grid origin </li>
<li>\ref get_ord : extract axis order array (\ref IPNT)</li>
</ul>
<li> Grid i/o package consists of functions
<ul>
<li>\ref read_grid : read grid from RSF header file </li>
<li>\ref par_grid : construct grid from \ref PARARRAY </li>
<li>\ref extend_array : extrapolate array data by constant along an axis </li>
<li>\ref rsfread : read \ref ireal array from RSF file pair </li>
<li>\ref rsfwrite : write \ref ireal array to RSF file pair </li>
</ul>
</ul>
<p>
Notes:
<p>
<ol>

<li>The implementation is independent of the Madagascar code package -
it uses only the RSF file format as described in the document
referenced above. 
</li>
<p>
<li>IWAVE::grid does <b>not</b> implement these two
features of Madagascar's package: 
<ul>
<li>DATAPATH: The output file (value associated to
keyword in) is a relative or fully qualified pathname in IWAVE::grid
usage, without any implicit completion. In particular, IWAVE::grid
<b>ignores</b> any DATAPATH environment variable, datapath=
specification on command line, datapath file in the working directory,
and so on.</li> 
<li>Only two data formats are currently recognized, namely native
4-byte floating point real and XDR-encoded 4-byte floating point real,
signalled by data_format=native_float and data_format=xdr_float
respectively. The esize keyword (a holdover from SEPlib) is
ignored. Thus all data is stored in 4-byte floating point form. The
data can however be read to, or written from, any format for which
implicit type conversion in C is defined, for example 8-byte real
(double). The the ireal typedef (see <a
href="../../../base/doc/html/index.html">IWAVE base library docs</a>
for more on this) defines the choice of internal data representation.
</li>
</ul>
</li> 
<p>
<li>IWAVE::grid provides two keyword choices in addition to those defined by Madagascar:
<ul>
<li>Value scaling: inclusion of the scale keyword causes the data to
be multiplied (divided) by the power of 10 indicated by the
corresponding value, as it is read (written). This feature permits the
user to change metric units 'on the fly', and avoid external
pre/postprocessing of data files. For example, to read data in m/ms
when it is stored on disk in m/s, insert 'scale=-3' in the RSF header
file.</li>
<li>Physical axis assignment: typical applications assign physical meaning to the coordinate axes (such as depth). To permit the axis ordering to be intperpreted properly without requiring that data files be transposed, IWAVE::grid provides keywords x_axis, y_axis, and z_axis to tag up to three axes in the data. The internal tags (axis id numbers) corresponding to these keywords are 2, 3, and 1 respectively. For example, including 'x_axis=1' in the RSF header file causes axis 1 (that is the fastest array index) to assigned with id=2. Applications using the grid data structure can use this information to identify axis 1 as 'the x axis', whatever that means in the application's terms. For example, for trace sampling it permits a sampling application to interpret shot and receiver x coordinates as indicating positions on axis 1.</li>
</ul>
</li>
<p>
<li>
Parallel implementation: all i/o takes place on rank
0. Implementation assumes that a chunk consisting of a number \ref
N_SLOWEST_AXIS_SLICE of (dim-1) slices may be safely stored in an
in-core buffer. Reads: data read into buffer from file on rank 0, then
broadcast to all other ranks. Each rank computes that part of the data
needed locally, then stores it. Writes: roughly the reverse of read,
but implemented with point-to-point communication, hence slower.
</li>
</ol>
<hr>
<a href="../../../doc/html/index.html">IWAVE Home Page</a>
*/
