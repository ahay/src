CC = gcc 
CFLAGS = -std=gnu9x -Wall -pedantic -I../lib -g
LDFLAGS = -L../lib -lrsfplot -lm

PLOTS = 7ab dc dcretard geqs geq2s headray vrms reflexpt reflkine

all: $(addsuffix .v, $(PLOTS))

%.x: %.c ../lib/librsfplot.a ../lib/rsfplot.h
	${CC} ${CFLAGS} $*.c ${LDFLAGS} -o $@

reflexpt.x: reflector.c ../lib/librsfplot.a ../lib/rsfplot.h
	${CC} ${CFLAGS} -DREFLEXPT reflector.c ${LDFLAGS} -o $@

reflkine.x: reflector.c ../lib/librsfplot.a ../lib/rsfplot.h
	${CC} ${CFLAGS} -DREFLKINE reflector.c ${LDFLAGS} -o $@

%.v: %.x
	$*.x > $@

geqs.v: sg.x
	sg.x 1. > $@

geq2s.v: sg.x
	sg.x 2. > $@

clean: 
	touch junk.x junk.v
	rm *.x *.v
