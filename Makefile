clean:
	cd seis/rsf; ${MAKE} distclean
	cd seis/main; ${MAKE} clean
	cd seis/proc; ${MAKE} clean
	cd seis/imag; ${MAKE} clean
	cd vplot/lib; ${MAKE} distclean
	cd vplot/main; ${MAKE} clean
