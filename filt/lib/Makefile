CC = gcc
CFLAGS = -Wall -pedantic -O3
CFLAGS = -Wall -pedantic -g -pg
CFLAGS = -std=gnu99 -Wall -pedantic -g
CFLAGS = -std=gnu9x -Wall -pedantic -g

DEPEND = -MM

SRC = alloc c99 decart error file files getpar math1 pfafft segy simtab stack
OBJS = $(addsuffix .o, $(SRC))	 
INCS = $(addsuffix .h, $(SRC)) 

all: librsf.a rsf.h

librsf.a: ${OBJS}
	ar rc $@ $^
	ranlib $@

rsf.h: ${INCS}
	cat ${INCS} > $@

%.o: %.c
	${CC} -c ${CFLAGS} $<

%.x: %.o ${OBJS}
	${CC} ${CFLAGS} $< ${OBJS} -lm -o $@

-include Makefile.dependencies

Makefile.dependencies:: *.c 
	$(CC) $(CFLAGS) $(DEPEND) $^ > Makefile.dependencies

clean:
	touch junk.o junk.x junk.H
	rm *.o *.x *.H

distclean:
	${MAKE} clean
	touch librsf.a rsf.h Makefile.dependencies
	rm librsf.a rsf.h Makefile.dependencies
