try:    from rsf.cluster import *
except: from rsf.proj    import *
    
# ------------------------------------------------------------
# DP test
def run(d,Lop,x,y,DEPEN=''):

    # x1
    Flow('_'+d+x+'1',x,'math output=0 | noise')
    # y1 = L  x1    
    Flow('_'+d+y+'1','_'+d+x+'1'+' '+DEPEN,Lop+' adj=n ')
    
    # y2
    Flow('_'+d+y+'2',y,'math output=0 | noise')
    # x2 = L' y2
    Flow('_'+d+x+'2','_'+d+y+'2'+' '+DEPEN,Lop+' adj=y ')

    for i in [x,y]:
        Flow('_'+d+i,['_'+d+i+'1','_'+d+i+'2'],
             '''
             add mode=p ${SOURCES[1]} |
             stack norm=n axis=1 | stack norm=n axis=1 | 
             stack norm=n axis=1 | stack norm=n axis=1
             ''')

    # DP values in percent
    Flow(d,['_'+d+x,'_'+d+y],
         'cat axis=1 ${SOURCES[1]} | scale axis=123 | scale rscale=1e2')
