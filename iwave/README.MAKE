[This file is deprecated, and will be removed in the next release. 
Please use instead the browser-enabled version 

iwave/doc/html/build.html 
]


           *** Using the TRIP Build System ***

0. Package structure

   The TRIP build system expects the files making up a package to be
   organized in a fashion common (though certainly not universal) in
   Un*x applications:

   include - a single directory (no subs), containing header files
   accessible across the package (and outside of it, as described in
   item 5 below). To function properly, headers should be named with
   one of the common extensions (.h for C headers, .hh, .H, .hpp for
   C++,...).

   src - source files, may have as many levels of subdirectory as you
   like, and contain header files used only locally. For proper
   compilation, use the common file extensions to designate these (.c
   for C, .cc, .C, or .cpp for C++, .f for F77, .f90 for f90,
   etc. etc.)
  
   lib - archive, tagged by package name and architecture tag (see
   1). If this directory is not present, no object archive will be
   created. If it is, all object files will be archived in
   lib/ARCH/lib<package>.a upon compilation via make. Thus: only one
   library per package, and static linking. Since the architecture tag
   is part of the path, the same package directory can contain
   libraries appropriate for several architectures, and the build
   system will link appropriately. [We may change this by-and-by so
   that the archive name depends on the compiler tag also - then archives
   built by different compiler families on the same architecture, as well
   as on different architectures, may be maintained in the same 
   directory tree. Not implemented as of 10.02.08.]

   main - if present, the install make target compiles and links any
   main program files (containing a line beginning "int main ...") in
   this directory. Note that the TRIP make rules use the suffix .x to
   designate executable files - whether made in main, or elsewhere. 
   [Commands defined by source in other directories are not made by
   the install target, but may be made by the default target. Only
   commands defined in main are "installed" in bin.]

   bin - if present, make will copy any command made in main to this
   directory, under an architecture-specific subdirectory.

   doc - not touched by the distribution version of the build system,
   common name for documentation directory.

1. Preparing the environment

   The TRIP build system combines make rules with C shell helper
   scripts. Since all of the scripts execute as csh subshells, you
   must define all environment variables parametrizing the build in a
   .cshrc in your home directory.

   The necessary rules and shells are contained in the TRIP admin package,
   which you should download from the TRIP subversion site, as described 
   on the TRIP software web page 

   http://www.trip.caam.rice.edu/software/software.html

   At least the following two variables must be defined in your .cshrc

   COMP: a tag identifying the compiler family
   TRIPADMIN: path to TRIP admin directory

   The compiler family tag is often, though not always, the name of the
   C compiler of the family. Its significance is as the last component 
   of the config file name (in TRIPADMIN/rules). 

   Run TRIPADMIN/shells/pvmgetarch to identify the architecture tag
   assigned to your platform, and examine (or create!) the config file
   TRIPADMIN/rules/Prg.defs.ARCH.COMP, where ARCH is the output of
   TRIPADMIN/shells/pvmgetarch. This file defines any non-default compiler
   and linker options. If any changes are needed, make them. Also make
   sure that your environment includes all information needed to
   configure the build. Some config files included with this release
   require environment variables to be defined beyond the two
   mentioned above. For example, our gcc config files parametrize the
   path to the compiler by GCCBIN, and our MPI-enabled config files
   identify the location of MPI via the variable MPIROOT. You will
   need to make sure that these or other environment parameters are
   defined in your .cshrc.

   For example, to build on OS X with the Gnu compilers, include 
 
   setenv COMP gcc
   setenv GCCBIN <path to gcc - often /usr/bin>
   
   in your .cshrc. TRIPADMIN/shells/pvmgetarch returns DARWIN when run on
   an OS X platform, so the config file included is

   TRIPADMIN/rules/Prg.defs.DARWIN.gcc.

   To build with MPI installed with the PGI suite on a 64 bit Linux 
   platform, include the lines

   setenv COMP mpi_pgcc
   setenv MPIROOT <path to the MPI installation you wish to use>

   in your .cshrc to use the rules in

   Prg.defs.LINUX64.mpi_pgcc

   Architecture/compiler choices available as this README is being edited:

   ARCH		COMP		description
   -------------------------------------------------------------------------------
   DARWIN	gcc		Mac OS-X, Gnu compiler suite
   DARWIN       mpi_gcc 	Max OS-X, MPI (mpich) compiled with gcc (mpicc etc.)
   LINUX	gcc		32 bit Linux, Gnu suite
   LINUX	mpi_gcc   	32 bit Linux, MPI (mpich) compiled with gcc (mpicc etc.)
   LINUX64	gcc		64 bit Linux, Gnu suite
   LINUX64	mpi_gcc   	64 bit Linux, MPI (mpich) compiled with gcc (mpicc etc)
   LINUX64      lam_gcc         64 bit Linux, MPI (lam) compiled with gcc (mpicc etc.)
   LINUX64	pgcc		64 bit Linux, PGI suite (pgcc etc.)
   LINUX64	mpi_pgcc   	64 bit Linux, MPI (mpich) compiled with pgcc (mpicc etc.)
   LINUX64      icc             64 bit Linux, Intel suite (icc etc.)
   LINUX64      mpi_icc         64 bit Linux, MPI (mpich) compiled with icc (mpicc etc.)
   LINUX64      mpi_altix       64 bit Linux, SGI MPI compiled with icc (mpicc etc.)

   Compiler Flags:

   The Prg.defs.ARCH.COMP file defines flags passed to the compiler
   and linker, usually as the value of the variable OCFLAGS. Three
   options are worth pointing out::

   -DSUXDR - flags inclusion of the XDR filters in various parts of the 
   IWAVE packages. [Chosen to toggle the XDR option in SU as well.]

   -DIWAVE_USE_MPI - set when MPI library calls should be included.

   -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -D_GNU_SOURCE - may be
   required to make gcc recognize large file (POSIX.1) extensions
   under some Linux distributions, eg. recent RedHats. Apparently
   defaulted under OS-X.

2. Instrumenting a package with makefiles

   In each package subdirectory, for example newwavem, if Makefile is not
   already in the directory, copy TRIPADMIN/makefile.boot to Makefile.
   
   Then creat makefiles throughout the package tree by

   make bootstrap

   Note that the install target (see below) depends on bootstrap if
   the makefile is the one copied from TRIPADMIN, but on the makefile
   target otherwise (after a package Makefile has been installed via
   bootstrap). After you have run bootstrap, the bootstrap target
   disappears. Normally you will not actually make bootstrap -
   instead, you will make install, which also compiles and archives
   object files and so on (see descriptions below)..

3. The make targets: 

   (default) - make all object and executable files in all subdirectories; 
   object files archived in the library lib/ARCH/lib<package>.a. Any
   executables compiled and (statically) linked to this library, and
   others as indicated by the hclasspath file (see below).

   install (available only in top level of package directory tree) -
   clean the directory tree (make clean, below); rebuild all makefiles
   throughout the tree (make makefile, or make bootstrap if top-level
   makefile lsstill in original state); compile all object files using
   headers in source directories, the package include directory, and
   any other header directories specified by the hclasspath (see
   below); archive object files in lib/ARCH/lib<package>.a, if a lib
   subdir is present; compile any commands in main (if a main subdir
   is present) and link these to lib/ARCH/lib<package>.a and to any
   external libraries; store copies of commands in bin subdir, if
   present.

   makefile - the makefile maintenance target. If you add or remove a
   file from somewhere in the package tree, execute this target at the
   top to rebuild makefiles throughout the tree appropriately. If the
   file in question is in a source directory with a local makefile,
   then make makefile in that directory is sufficient.

   lib - make all object files, archive in lib/ARCH/lib<package>.a

   clean - remove object and executable files throughout tree (except
   for executable copies in bin, see below), other secondary files
   as defined in clean rules (TRIPADMIN/rules/Doc.defs.clean).
 
   libclean - clean plus remove length-0 files used for their timestamps
   to keep libs up-to-date. Should use this target if you want to force
   regeneration of the library.
 
   burn - clean, then remove object archive (in lib/ARCH) and all makefiles.
   This target returns the package to its original, pristine condition.

   Besides the native suffixes (.c etc. for source files, .o for object
   files), 

4. Executable maintenance

   If main and bin subdirectories are present, all executables made
   in main are copied to bin/ARCH. The install target makes all of these 
   and copies them to bin/ARCH. The copies in bin/ARCH are not deleted 
   by clean, but are deleted by burn. Thus make clean cleans up the
   directory but leaves the command created from it available (in
   bin); make burn gets rid of the commands as well.

5. Linking to external libraries

   To link libraries and headers from outside the subpackage tree,
   create a file named hclasspath in the root directory of the
   subpackage, and write into it a path, using exactly the syntax used
   in defining environment paths: for example

   $FOO:../foobar

   links to libraries, and uses headers, in FOO and ../foobar (FOO being
   defined in the environment). The list must respect dependencies, 
   exactly as in a linker command: if FOO depends on ../foobar, then 
   the two must be listed in the order shown above. To link properly, 
   an external package must be structured in the same way as those 
   built by the TRIP tools, that is: each package, say FOO, must 
   define precisely one library (archive), stored (for example) in
   $FOO/lib/$ARCH/lib$FOO.a (note the architecture subdir). Header
   files (to be made accessible to the compiler) must reside in an
   include subdirectory under the root of the foreign package, for
   example $FOO/include. If the library you wish to use does not have 
   this structure to begin with, you may need to create a link tree 
   mimicing it.

6. Credit
 
   The TRIP Build System is based on an early version of the Stanford
   Exploration Project make rules. The makefile autowriter maw and
   other automation scripts were developed by TRIP. None of this would
   be possible without the Free Software Foundation and its Gnu make
   utility, several features of which are critical to the success of
   our build scheme.

WWS, 01.08
