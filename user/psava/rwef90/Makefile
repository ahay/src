include ${SEPINC}/SEP.top
MAKE_DEPEND=yes
SAVE_FORTRAN=no
include ${PCSINC}/PCS.prg

SRCDIR=./src
BINDIR=.
LIBDIR=.
OBJDIR=.
INCDIR=.

UF90INCLUDES=-I.
OLEVEL=2
UF90LIBS = -Bstatic
UF90FLAGS= -O2 #-Mbounds

default: $B/RWEmig.x $B/RWEab.x $B/RWEaj.x $B/CC2RC.x

clean:jclean
include ${SEPINC}/SEP.bottom

