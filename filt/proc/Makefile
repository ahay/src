CFLAGS = -O3 
CFLAGS = -std=gnu99 -Wall -g -pedantic -I../rsf
LDFLAGS = -L../rsf -lrsf -lm
DEPEND = -MM

PROGS = bandpass bin bin1 depth2time dip envelope nmo noise \
	pwd randrefl sigmoid smooth spectra taup time2depth warp1

bandpass.OBJS = butterworth
bin.OBJS = int2 interp adjnull
bin1.OBJS = int1 interp adjnull
depth2time.OBJS = stretch tridiagonal
dip.OBJS = dip3 divide1 allp3 tridiagonal apfilt
envelope.OBJS = triangle
nmo.OBJS = stretch tridiagonal
noise.OBJS = randn
pwd.OBJS = allp3 apfilt
randrefl.OBJS = randn
sigmoid.OBJS = random
smooth.OBJS = triangle
spectra.OBJS = empty
taup.OBJS = radon cadj ctoeplitz
time2depth.OBJS = fint1 tridiagonal spline interp_spline extend banded
warp1.OBJS = divlap1 int1 interp_spline prefilter recfilt adjnull banded

MAINS = $(addprefix sf, $(PROGS))
OBJS = $(addsuffix .o, ${$*.OBJS})

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
