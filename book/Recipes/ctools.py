try:    from rsf.cluster import *
except: from rsf.proj    import *

# ------------------------------------------------------------
def Temp(o,i,r):
    Flow(o,i,r+ ' datapath=%s '%os.environ.get('TMPDATAPATH'))
# ------------------------------------------------------------

# ------------------------------------------------------------
# prepare cluster run
def cls(job,    # job name
        nodes): # number of cluster nodes

    cluster=dict()

    cluster['job']=job
    cluster['user']=os.environ.get('USER')
    cluster['email']=cluster['user']+'@mines.edu'
    cluster['time']=5
    cluster['ppn']=12
    cluster['nodes']=nodes

    return cluster
    
# ------------------------------------------------------------
# prepare groups for cluster run
def grp(nf,of,df, # input files n,o,d
        ng):      # number of files in a group

    groups=[]
    for im in range(int(nf/ng)+1): # loop over groups

        mlo=    of+(im  )*ng*df           # LOw  group index
        mhi=min(of+(im+1)*ng*df,of+nf*df) # HIgh group index
        
        if(mlo != mhi):
            groups.append(range(mlo,mhi,df)) # group indices

    return(groups)

# ------------------------------------------------------------
# cat "nf" files in groups of "ng" files
def cat(output,   # output
        prefix,   # input template
        nf,of,df, # input files n,o,d
        ng,       # number of files in a group
        axis,     # cat axis
        time=5,   # cluster execution time
        nodes=1): # number of cluster nodes

    # sort groups
    groups=grp(nf,of,df,ng)
    nodes=min(nodes,len(groups))

    # loop over groups
    Fork(time=time,ipn=len(groups)/nodes,nodes=nodes)
    for ig in range(len(groups)):
        Flow(output+"%04d"%ig,
             [prefix%jg for jg in groups[ig]],
             'cat axis=%d space=n ${SOURCES[1:%d]}'%(axis,len(groups[ig])))
        Iterate()
    Join()

    # join groups
    Flow(output,
         [output+"%04d"%ig for ig in range(len(groups))],
         'cat axis=%d space=n ${SOURCES[1:%d]}'%(axis,len(groups)))

# ------------------------------------------------------------
# add "nf" files in groups of "ng" files
def add(output,   # output
        prefix,   # input template
        nf,of,df, # input files n,o,d
        ng,       # number of files in a group
        time=5,   # cluster execution time
        nodes=1): # number of cluster nodes

    # sort groups
    groups=grp(nf,of,df,ng)
    nodes=min(nodes,len(groups))

    # loop over groups
    Fork(time=time,ipn=len(groups)/nodes,nodes=nodes)
    for ig in range(len(groups)):
        Flow(output+"%04d"%ig,
             [prefix%jg for jg in groups[ig]],
             'add ${SOURCES[1:%d]}'%(len(groups[ig])))
        Iterate()
    Join()
        
    # join groups
    Flow(output,
         [output+"%04d"%ig for ig in range(len(groups))],
         'add ${SOURCES[1:%d]}'%(len(groups)))
