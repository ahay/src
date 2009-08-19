from pylab import *
from numpy import *

down = '''
54
166
279
357
396
469
531
582
646
789
927
1077
1187
1281
1368
1463
1662
1837
1929
2038
2109
2236
2352
2539
2705
2863
2993
3146
3340
3731
4072
4344
4668
5070
5258
5452
5765
'''
down = map(int,down.split())
ind = arange(len(down))
width=0.75

p1 = bar(ind,down,width)
ylabel('Downloads')
xticks(ind[::6]+width/2, ('Jun-2006', 'Dec-2006',
                          'Jun-2007', 'Dec-2007', 
                          'Jun-2008', 'Dec-2008'))
title('Cumulative Downloads')
#grid(True)

#show()
savefig('junk_py.eps')



