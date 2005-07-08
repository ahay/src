#!/bin/env python

import os,commands, re

retime = re.compile('(\w+)\s([\d\.]+)m([\d\.]+)s')
out = open('fft.rsf','w')

for size in xrange(1,50001):
    spike = 'spike%d.rsf' % size
    fft = 'fft%d.rsf' % size
    os.system('sfspike n1=%d n2=1000 > %s' % (size,spike))
    times = commands.getoutput('time sffft1 < %s > %s' % (spike,fft))
    total = 0.0
    for time in times.split('\n'):
        got = retime.match(time)
        if got and (got.group(1)=='user' or got.group(1)=='sys'):
            total += 60*float(got.group(2))+float(got.group(3))
    print 'time[%d]=%g' % (size,total)
    out.write('%g\n' % total)
    os.system('sfrm %s %s' % (spike,fft))

out.write('n1=%d o1=1 d1=1\n' % size)
out.write('in=fft.rsf\n')
out.write('data_format=ascii_float\n')
out.close()

