TARGET = petsc_test
include $(PETSC_DIR)/conf/variables

default: test
all: default

test.o: _test.c 
	${PCC} ${PCC_FLAGS} ${CCPPFLAGS} -c $< -o $@

test: test.o
	$(CLINKER) -o $@ $< $(PETSC_LIB)

clean:
	-rm -f test test.o 

.PRECIOUS: test test.o
