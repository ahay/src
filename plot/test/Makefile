CC = gcc 
CFLAGS = -std=gnu9x -Wall -pedantic -I../lib -g
LDFLAGS = -L../lib -lrsfplot -lm

%.x: %.c ../lib/librsfplot.a ../lib/rsfplot.h
	${CC} ${CFLAGS} $*.c ${LDFLAGS} -o $@

geqs.v: sg.x
	sg.x 1. > $@

geq2s.v: sg.x
	sg.x 2. > $@

headray.v: headray.x
	headray.x > $@

vrms.v: vrms.x
	vrms.x > $@

clean: 
	touch junk.x junk.v
	rm *.x *.v
