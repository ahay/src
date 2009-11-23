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
6021
6303
6665
7038
'''
down = map(int,down.split())
ind = arange(len(down))
width=0.75

p1 = bar(ind,down,width)
ylabel('Downloads')
xticks(ind[::6]+width/2, ('06/2006', '12/2006',
                          '06/2007', '12/2007', 
                          '06/2008', '12/2008',
                          '06/2009'))
title('Cumulative Downloads')
#grid(True)

#show()
savefig('junk_py.eps')



