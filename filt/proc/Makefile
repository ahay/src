CFLAGS = -O3 
CFLAGS = -Wall -g -pedantic -I../rsf
LDFLAGS = -L../rsf -lrsf -lm
DEPEND = -MM

PROGS = bin bin1 dip pwd sigmoid smooth

bin.OBJS = int2 interp adjnull
bin1.OBJS = int1 interp adjnull
dip.OBJS = dip3 divide1 allp3 tridiagonal apfilt
pwd.OBJS = allp3 apfilt
sigmoid.OBJS = random
smooth.OBJS = triangle

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
	touch sfjunk junk.H junk.lint junk.o junklint 
	rm sf* *.H *.lint *.o junklint
