CFLAGS = -O3 
CFLAGS = -std=gnu99 -Wall -g -pedantic -I../rsf
LDFLAGS = -L../rsf -lrsf -lm
DEPEND = -MM

PROGS = cell2 eikonal gazdag rays2 sag shoot2 symes voft

cell2.OBJS = celltrace cell eno2 eno quadratic
eikonal.OBJS = fastmarch pqueue neighbors
gazdag.OBJS = gazdag
rays2.OBJS = raytracing atela grid3 grid2 eno3 eno2 eno
sag.OBJS = empty 
shoot2.OBJS = celltrace cell eno2 eno quadratic fzero
symes.OBJS = symes
voft.OBJS = empty

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
