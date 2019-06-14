SConstruct.jupyter: *.scons
	echo "from rsf.proj import *" > $@
	cat $^ >> $@
	echo "" >> $@
	echo "End()" >> $@

%.png: SConstruct.jupyter
	scons -f SConstruct.jupyter Fig/$*.vpl
	vpconvert pen=gd fat=3 serifs=n bgcolor=w Fig/$*.vpl $@
