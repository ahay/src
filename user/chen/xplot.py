'extra plots'

from rsf.proj import *

def Mplot(target, source):
	Result(target, source,
		'''
		pspen serifs=y fat=3 color=y label="" tex=y scale=0.8
		| ps2pdf - -
		| pdfjam -q -o /dev/stdout --papersize "{5in,5in}"
		--trim '6.2cm 2.2cm 7.2cm 7.5cm'
		| pdfjam -q -o /dev/stdout --papersize "{5in,5in}"
		--nup 3x3 
		''',
		suffix='.pdf', src_suffix='.vpl')


