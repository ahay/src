CFLAGS = -O3 
CFLAGS = -std=gnu9x -Wall -g -pedantic -I../rsf
LDFLAGS = -L../rsf -lrsf -lm
DEPEND = -MM

PROGS = bandpass bin bin1 depth2time dip envelope fft1 helicon invbin1 nmo \
	noise pwd randrefl reshape ricker sigmoid smooth spectra taup taupmo \
	time2depth warp1 window1

bandpass = butterworth
bin = int2 interp adjnull
bin1 = int1 interp adjnull
depth2time = stretch tridiagonal
dip = dip3 divide1 allp3 tridiagonal apfilt
envelope = triangle
fft1 = empty
helicon = polydiv regrid helicon helix adjnull
invbin1 = int1 interp adjnull chain bigsolver tcai1 cgstep causint \
	triangle1 triangle
nmo = fint1 tridiagonal spline interp_spline extend banded
noise = randn
pwd = allp3 apfilt
randrefl = randn
reshape = empty
ricker = empty
sigmoid = random
smooth = triangle
spectra = empty
taup = radon cadj ctoeplitz
taupmo = fint1 tridiagonal spline interp_spline extend banded
time2depth = fint1 tridiagonal spline interp_spline extend banded
warp1 = divlap1 int1 interp_spline prefilter recfilt adjnull banded
window1 = window1

MAINS = $(addprefix sf, $(PROGS))
OBJS = $(addsuffix .o, ${$*})

all: ${MAINS}

%.o: %.c
	${CC} -c ${CFLAGS} $<

sf%: M%.o ../rsf/librsf.a ../rsf/rsf.h 
	${MAKE} ${OBJS}
	${CC} ${CFLAGS} M$*.o ${OBJS} ${LDFLAGS} -o $@

%.lint: %.c ../rsf/rsf.h
	-splint +posixlib -I../rsf $^ > $@

junklint: $(addsuffix .lint, $(PROGS))
	cat $^ > $@

-include Makefile.dependencies

Makefile.dependencies:: *.c 
	$(CC) $(CFLAGS) $(DEPEND) $^ > Makefile.dependencies

clean:
	touch sfjunk junk.H junk.lint junk.o junklint Makefile.dependencies
	rm sf* *.H *.lint *.o junklint Makefile.dependencies
