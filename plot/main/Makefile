CC = cc
CC = gcc 

RSF = ../../seis/rsf
RSFPLOT = ../lib

CFLAGS = -I${RSF} -I${RSFPLOT} -g
CFLAGS = -std=gnu9x -Wall -pedantic -I${RSF} -I${RSFPLOT} -g

LDFLAGS =  -L${RSFPLOT} -lrsfplot -L${RSF} -lrsf -lm -lnsl
LDFLAGS =  -L${RSFPLOT} -lrsfplot -L${RSF} -lrsf -lm

PROGS = contour dots graph grey plas pldb plotrays

MAINS = $(addprefix sf, $(PROGS))

all: ${MAINS}

sf%: %.c ${RSF}/rsf.h ${RSF}/librsf.a \
	${RSFPLOT}/rsfplot.h ${RSFPLOT}/librsfplot.a
	${CC} ${CFLAGS} $*.c ${LDFLAGS} -o $@

clean:
	touch sfjunk junk.H junk.v junk.lint junklint 
	rm sf* *.H *.v *.lint junklint
