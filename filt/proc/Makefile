CC = cc
CC = gcc
CFLAGS = -O3 
CFLAGS = -g -I../rsf
CFLAGS = -std=gnu9x -Wall -g -pedantic -I../rsf

LDFLAGS = -L../rsf -lrsf -lm
LDFLAGS = -L../rsf -lrsf -lm -lnsl

DEPEND = -xM1
DEPEND = -MM

PROGS = ai2refl bandpass bin bin1 depth2time dip envelope fft1 fft2 fft3 \
	helicon invbin1 invwarp nmo noise pef pwd randrefl reshape ricker \
	sigmoid smooth spectra stretch2 taup taupmo time2depth warp1 \
	window1 window2 xcorr1

ai2refl = empty
bandpass = butterworth
bin = int2 interp adjnull
bin1 = int1 interp adjnull
depth2time = stretch tridiagonal
dip = dip3 divide1 allp3 tridiagonal apfilt
envelope = triangle
fft1 = empty
fft2 = empty
fft3 = empty
helicon = polydiv regrid helicon helix adjnull
invbin1 = int1 interp adjnull chain bigsolver tcai1 cgstep causint \
	triangle1 triangle polydiv pef bound helix hconest helicon regrid
invwarp = stretch tridiagonal
nmo = fint1 tridiagonal spline interp_spline extend banded
noise = randn
pef = pef printfilter boxfilter helix bound misinput regrid chain \
	createhelix compress helicon hconest cgstep bigsolver adjnull
pwd = allp3 apfilt
randrefl = randn
reshape = empty
ricker = empty
sigmoid = random
smooth = triangle
spectra = empty
stretch2 = stretch2 banded
taup = radon cadj ctoeplitz
taupmo = fint1 tridiagonal spline interp_spline extend banded
time2depth = fint1 tridiagonal spline interp_spline extend banded
warp1 = divlap1 int1 interp_spline prefilter recfilt adjnull banded
window1 = window1
window2 = window2
xcorr1 = window1 xcorr stretch2 banded

MAINS = $(addprefix sf, $(PROGS))
OBJS = $(addsuffix .o, ${$*})

all: ${MAINS}

%.o: %.c
	${CC} -c ${CFLAGS} $<

sf%: M%.o ../rsf/librsf.a ../rsf/rsf.h 
	${MAKE} ${OBJS}
	${CC} ${CFLAGS} M$*.o ${OBJS} ${LDFLAGS} -o $@

Test%.x: Test%.c %.o
	${CC} ${CFLAGS} Test$*.c $*.o ${LDFLAGS} -o $@ 

%.lint: %.c ../rsf/rsf.h
	-splint +posixlib -I../rsf $^ > $@

junklint: $(addsuffix .lint, $(PROGS))
	cat $^ > $@

-include Makefile.dependencies

Makefile.dependencies:: *.c 
	$(CC) $(CFLAGS) $(DEPEND) $^ > Makefile.dependencies

clean:
	touch sfjunk junk.H junk.lint junk.o junklint Makefile.dependencies Testjunk.x
	rm sf* *.H *.lint *.o junklint Makefile.dependencies Test*.x
