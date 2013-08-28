################################################################################
###
### Copyright Rice University, 2004.
### All rights reserved.
###
### Permission is hereby granted, free of charge, to any person obtaining a
### copy of this software and associated documentation files (the "Software"),
### to deal in the Software without restriction, including without limitation
### the rights to use, copy, modify, merge, publish, distribute, and/or sell
### copies of the Software, and to permit persons to whom the Software is
### furnished to do so, provided that the above copyright notice(s) and this
### permission notice appear in all copies of the Software and that both the
### above copyright notice(s) and this permission notice appear in supporting
### documentation.
### 
### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
### IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
### FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
### RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
### NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
### DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
### PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
### ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
### THIS SOFTWARE.
###
### Except as contained in this notice, the name of a copyright holder shall
### not be used in advertising or otherwise to promote the sale, use or other
### dealings in this Software without prior written authorization of the
### copyright holder.
###
###############################################################################
 
########################## LOCAL RULES ############################
 
PACKAGE    = demo_new
SHELL      = /bin/sh
TOP        = ./..
HERE       = ./
RULES      = ${TRIPADMIN}/rules
SHELLS     = ${TRIPADMIN}/shells
LASTARCH   = LINUX64
LASTENDIAN = little
TRIPARCH   = $(shell (exec ${SHELLS}/pvmgetarch))
ENDIAN     = $(shell (exec ${SHELLS}/getendian))
SRCROOT    = .
DATAROOT   = .
BINDIR     =
INCLUDES   = ${IWAVEPP}/asg++/include ${IWAVEPP}/iwave++/include ${IWAVEPP}/grid++/include ${IWAVEPP}/segy++/include ${IWAVEPP}/fos/include ${IWAVEPP}/oc/include ${RVL}/tsopt/include ${RVL}/umin/include ${RVL}/alg/include ${RVL}/mpirvl/include ${RVL}/local/include ${RVL}/rvl/include ${RVL}/rand/include ${IWAVE}/asg/include ${IWAVE}/sample/include ${IWAVE}/iwave/include ${IWAVE}/grid/include ${IWAVE}/trace/include ${IWAVE}/sucore/include ${IWAVE}/base/include
 
LIBDIR =
 
############################# TOP STANDARD RULES #######################
 
include ${RULES}/Make.top
 
############################# LOCAL DEFINITIONS ########################
 
LIB = 
 
LNKLIB = ${IWAVEPP}/asg++/lib/${TRIPARCH}/libasg++.a ${IWAVEPP}/iwave++/lib/${TRIPARCH}/libiwave++.a ${IWAVEPP}/segy++/lib/${TRIPARCH}/libsegy++.a ${RVL}/tsopt/lib/${TRIPARCH}/libtsopt.a ${RVL}/umin/lib/${TRIPARCH}/libumin.a ${RVL}/mpirvl/lib/${TRIPARCH}/libmpirvl.a ${RVL}/rand/lib/${TRIPARCH}/librand.a ${IWAVE}/asg/lib/${TRIPARCH}/libasg.a ${IWAVE}/sample/lib/${TRIPARCH}/libsample.a ${IWAVE}/iwave/lib/${TRIPARCH}/libiwave.a ${IWAVE}/grid/lib/${TRIPARCH}/libgrid.a ${IWAVE}/trace/lib/${TRIPARCH}/libtrace.a ${IWAVE}/sucore/lib/${TRIPARCH}/libsucore.a ${IWAVE}/base/lib/${TRIPARCH}/libbase.a 
UCLIB   = ${LNKLIB}
UC++LIB = ${LNKLIB}
UF77LIB = ${LNKLIB}
UF90LIB = ${LNKLIB}
 
 
RESULTSER = rami
 
######################## TARGET DEFINITIONS ############################
 
default:
 
 
burn: distclean
	@(echo "barf" > /dev/null)
 
clean: jclean
 
libclean:
	@(echo "barf" > /dev/null)
 
makefile: makefile.local makefile.subs.force
 
makefile.local:
	@(${SHELLS}/maw ${TOP})
 
makefile.subs.force:
	@(${MAKE} makefile.subs -f ./Makefile)
 
makefile.subs:
	@(echo "barf" > /dev/null)
 
checkarch:
	@(if [ ${TRIPARCH} != ${LASTARCH} ]; then \
		(${MAKE} makefile.local; \
		${MAKE} jclean); \
	fi )
 
include tests.er
 
######################### BOTTOM STANDARD RULES ########################
 
include $(RULES)/Make.bottom
 
####################### HEADER DEPENDENCE RULES ########################
 
