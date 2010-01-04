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
7392
7824
'''
down = map(int,down.split())
ind = arange(len(down))
width=0.75

p1 = bar(ind,down,width)
ylabel('Downloads')
xticks(ind[::6]+width/2, ('06/06', '12/06',
                          '06/07', '12/07', 
                          '06/08', '12/08',
                          '06/09', '12/09'))
title('Cumulative Downloads')
#grid(True)

#show()
savefig('junk_py.eps')



