include ${PSP_DIR}/conf/pspvariables

options:
	@echo ${CXX} ${PSP_COMPILE_FLAGS} ${PSP_LINK_FLAGS} ${PSP_LIBS}

# Since not all libs are linked with "-l", this avoids the problem of parsing
# the libraries from the above command
libs:
	@echo ${PSP_LIBS}
