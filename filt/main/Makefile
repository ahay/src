CC = gcc 
CFLAGS = -Wall -pedantic -I../rsf -O3
CFLAGS = -Wall -pedantic -I../rsf -g -pg
CFLAGS = -std=gnu99 -Wall -pedantic -I../rsf -g
CFLAGS = -std=gnu9x -Wall -pedantic -I../rsf -g

LDFLAGS = -L../rsf -lrsf -lm

PROGS = add attr cat cmplx cp dd disfil get headermath imag in interleave \
	log math merge mv pad real reverse rm rtoc scale \
	segyread segywrite spike stack transp window

MAINS = $(addprefix sf, $(PROGS))

all: ${MAINS}

sf%: %.c ../rsf/librsf.a ../rsf/rsf.h
	${CC} ${CFLAGS} $*.c ${LDFLAGS} -o $@

%.lint: %.c ../rsf/rsf.h
	-splint +posixlib -I../rsf $^ > $@

junklint: $(addsuffix .lint, $(PROGS))
	cat $^ > $@

clean:
	touch sfjunk junk.H junk.lint junklint 
	rm sf* *.H *.lint junklint
