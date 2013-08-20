/**
@page parallel Parallel Simulations with IWAVE 
<p> 

IWAVE simulations can run in either serial (single thread) mode, or
distributed (MPI) or shared-memory (OpenMP) parallel modes. Serial
mode is the default. For more on installation of the parallel options,
see the IWAVE <a href="../../../README.INSTALL">install </a> and <a
href="../../../README.MAKE">build </a> notes.

<p> 
Distributed mode
supports either loop-level (domain decomposition) or task level
(parallel shots) parallelism, or both simultaneously.  
<p> 

<b>Cartesian domain decomposition via MPI</b>: defined by specifying the number of
domains along each coordinate axis. 
<ul> 
<li> mpi_np1 = [int] [Default = 1] number of subdomains along axis 1</li>
<li> mpi_np2 = [int] [Default = 1] number of subdomains along axis 2</li>
<li> mpi_np3 = [int] [Default = 1].number of subdomains along axis 3</li>
</ul>

In domain decomposition mode, IWAVE divides regular rectangular
simulation grids into a rectangular array of regular regular
rectangular subgrids, of as nearly equal dimensions as possible. Thus
mpi_np2=2 indicates that a 100 x 100 x 100 grid will be partitioned
into two subgrids, each of size 100 x 50 x 100 (roughly - the actual
number may be influenced by added ghost cells needed for data exchange
at subdomain boundaries, as determined by scheme stenceils, or by PML
or other boundary conditions). Setting mpi_np1=2, mpi_np2=2 will yield
four domains of approximate size 50 x 50 x 100, and so on. The code
currently makes no attempt to account for load imbalance, for
instance, between PML and non-PML subgrids. The number of subdomains
can only be set > 1 in any direction if MPI parallelism is installed,
and if the spatial dimension of the problem is at least the index of
that direction - that is, the code exits with an error message if
mpi_np3 > 1 and the problem spatial dimnension is 1 or 2. The number
of MPI processes specified in the command line or batch script must be
at least the total number of domains = mpi_np1*mpi_np2*mpi_np3.  <p>

<b>Task parallelism</b>, that is, simultaneous simulation of several
shots: toggled by a nonzero value of the partask parameter:

<ul> 
<li> partask = [int] [Default = 0] if set, load as many
as possible given the size of MPI_COMM_WORLD,
<i> but at most partask</i>, simulations.</li> 
</ul>

For example, if mpi_np1=mpi_np2=mpi_np3=2 and 32 processors are assigned to the simulation, then setting partask to be at least 4 will cause up to 4 shots to be loaded at a time, simulated, and written to the output file, until no shots remain to be simulated. If partask is less than 4, but greater than zero, then partask simulations will be loaded. If 33 processors were assigned and partask is at least 4, IWAVE will still load 4 shots at a time, but the 33rd process will be flagged as inactive, that is, not used in the simulation. If 31 processors are assigned and partask is at least three, then IWAVE loads 3 shots at a time and flags the 25th-31st processes as inactive.
<p>
It is possible to use domain decomposition without task parallelism (partask=0 or number of processes less than twice the number of domains), or task parallelism without domain decomposition (partask > 0, mpi_np1=mpi_np2=mpi_np3=1), or both simultaneously.
<p>
Note that the user must define an upper bound (partask) on the number of simultaneous simulations. The actual number will be the minimum of this user-defined upper bound and the maximum that can fit in MPI_WORLD_COMM. This user role is unavoidable, as the assignment of processes takes place before the simulation geometry data is read.
<p>
<b>Threaded parallelism via OpenMP</b>: Regulated by 
<ul>
<li>omp_nt = [int] [Default = 1] number of OpenMP threads -
significant only if OpenMP is enabled.during install</li> 
</ul> 
<p>
Use requires that stencil code be instrumented with appropriate pragmas. Have kept this option available for possible future use. With current CPU designs, multithreaded execution of standard finite difference stencils does not seem to result in useful speedup.
*/
