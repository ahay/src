CC = gcc
CFLAGS = -Wall -pedantic -O3
CFLAGS = -Wall -pedantic -g -pg
CFLAGS = -std=gnu99 -Wall -pedantic -g
CFLAGS = -std=gnu9x -Wall -pedantic -g -I../../seis/rsf
LDFLAGS = -L../../seis/rsf -lrsf

DEPEND = -MM

SRC = axis coltab gplot stdplot vplot
OBJS = $(addsuffix .o, $(SRC))	 
INCS = $(addsuffix .h, $(SRC)) 

all: librsfplot.a rsfplot.h

librsfplot.a: ${OBJS}
	ar rc $@ $^
	ranlib $@

rsfplot.h: ${INCS}
	cat ${INCS} > $@

%.o: %.c
	${CC} -c ${CFLAGS} $<

%.x: %.o ${OBJS}
	${CC} ${CFLAGS} $< ${OBJS} ${LDFLAGS} -lm -o $@ 

-include Makefile.dependencies

Makefile.dependencies:: *.c 
	$(CC) $(CFLAGS) $(DEPEND) $^ > Makefile.dependencies

clean:
	touch junk.o junk.x junk.H
	rm *.o *.x *.H

distclean:
	${MAKE} clean
	touch librsfplot.a rsfplot.h Makefile.dependencies
	rm librsfplot.a rsfplot.h Makefile.dependencies
