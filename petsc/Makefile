TARGET = petsc_test
include $(PETSC_DIR)/conf/variables

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.c, %.o, $(wildcard *.c))
HEADERS = $(wildcard *.h)

%.o: %.c $(HEADERS)
	${PCC} ${PCC_FLAGS} ${CCPPFLAGS} -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CLINKER) -o $@ $(OBJECTS) $(PETSC_LIB)

clean:
	-rm -f *.o
	-rm -f $(TARGET)
