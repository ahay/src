CFLAGS = -O3 
CFLAGS = -std=gnu99 -Wall -g -pedantic -I../rsf
LDFLAGS = -L../rsf -lrsf -lm
DEPEND = -MM

PROGS = cell2 eikonal enoder1 gazdag gauss pick0 pickd \
	rays2 sag shoot2 sstep1 symes voft

cell2 = celltrace cell eno2 eno quadratic
eikonal = fastmarch pqueue neighbors
enoder1 = eno
gazdag = gazdag
gauss = empty
pick0 = pick0 eno
pickd = pick0 eno
rays2 = raytracing atela grid3 grid2 eno3 eno2 eno
sag = empty 
shoot2 = celltrace cell eno2 eno quadratic fzero
sstep1 = split1
symes = symes
voft = empty

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
