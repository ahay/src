SConstruct: *.scons
	echo "from rsf.proj import *" > $@
	cat $^ >> $@
	echo "" >> $@
	echo "End()" >> $@

%.png: SConstruct
	scons Fig/$*.vpl
	vpconvert pen=gd fat=3 serifs=n bgcolor=w Fig/$*.vpl $@
