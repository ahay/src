try:    from rsf.cluster import *
except: from rsf.proj    import *

# ------------------------------------------------------------
def Temp(o,i,r):
    Flow(o,i,r+ ' datapath=%s '%os.environ.get('TMPDATAPATH',os.environ.get('DATAPATH')))
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
def grplist(list, # list of input files 
            ng):  # number of files in a group

    nf=len(list)
    of=0
    df=1

    groups=[]
    for im in range(int(nf/ng)+1): # loop over groups

        mlo=    of+(im  )*ng*df           # LOw  group index
        mhi=min(of+(im+1)*ng*df,of+nf*df) # HIgh group index
        
        if(mlo != mhi):
            groups.append(list[mlo:mhi]) # group indices

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
    if(nodes>1): Fork(time=time,ipn=len(groups)/nodes,nodes=nodes)
    for ig in range(len(groups)):
        Flow(output+"-g%04d"%ig,
             [prefix%jg for jg in groups[ig]],
             'cat axis=%d space=n ${SOURCES[1:%d]}'%(axis,len(groups[ig])))
        
        if(nodes>1): Iterate()
    if(nodes>1): Join()

    # join groups
    Flow(output,
         [output+"-g%04d"%ig for ig in range(len(groups))],
         'cat axis=%d space=n ${SOURCES[1:%d]}'%(axis,len(groups)))

def catlist(output,   # output
            prefix,   # input template
            list,
            ng,       # number of files in a group
            axis,     # cat axis
            time=5,   # cluster execution time
            nodes=1): # number of cluster nodes

    # sort groups
    groups=grplist(list,ng)
    nodes=min(nodes,len(groups))
    
    # loop over groups
    if(nodes>1): Fork(time=time,ipn=len(groups)/nodes,nodes=nodes)
    for ig in range(len(groups)):
        Flow(output+"-g%04d"%ig,
             [prefix%jg for jg in groups[ig]],
             'cat axis=%d space=n ${SOURCES[1:%d]}'%(axis,len(groups[ig])))
        
        if(nodes>1): Iterate()
    if(nodes>1): Join()

    # join groups
    Flow(output,
         [output+"-g%04d"%ig for ig in range(len(groups))],
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
    if(nodes>1): Fork(time=time,ipn=len(groups)/nodes,nodes=nodes)
    for ig in range(len(groups)):
        Flow(output+"-g%04d"%ig,
             [prefix%jg for jg in groups[ig]],
             'add ${SOURCES[1:%d]}'%(len(groups[ig])))
        if(nodes>1): Iterate()
    if(nodes>1): Join()
        
    # join groups
    Flow(output,
         [output+"-g%04d"%ig for ig in range(len(groups))],
         'add ${SOURCES[1:%d]}'%(len(groups)))

def addlist(output,   # output
            prefix,   # input template
            list,
            ng,       # number of files in a group
            time=5,   # cluster execution time
            nodes=1): # number of cluster nodes

    # sort groups
    groups=grplist(list,ng)
    nodes=min(nodes,len(groups))

    # loop over groups
    if(nodes>1): Fork(time=time,ipn=len(groups)/nodes,nodes=nodes)
    for ig in range(len(groups)):
        Flow(output+"-g%04d"%ig,
             [prefix%jg for jg in groups[ig]],
             'add ${SOURCES[1:%d]}'%(len(groups[ig])))
        if(nodes>1): Iterate()
    if(nodes>1): Join()
        
    # join groups
    Flow(output,
         [output+"-g%04d"%ig for ig in range(len(groups))],
         'add ${SOURCES[1:%d]}'%(len(groups)))


# ------------------------------------------------------------
# make movie of "nf" plots in groups of "ng" files
def mov(output,   # output
        prefix,   # input template
        nf,of,df, # input files n,o,d
        ng,       # number of files in a group
        time=5,   # cluster execution time
        nodes=1): # number of cluster nodes

    # sort groups
    groups=grp(nf,of,df,ng)
    nodes=min(nodes,len(groups))

    # loop over groups
    if(nodes>1): Fork(time=time,ipn=len(groups)/nodes,nodes=nodes)
    for ig in range(len(groups)):
        Plot(output+"-g%04d"%ig,
             [prefix%jg for jg in groups[ig]],
             'Movie')
        if(nodes>1): Iterate()
    if(nodes>1): Join()

    # join groups
    Result(output,
         [output+"-g%04d"%ig for ig in range(len(groups))],
         'Movie')

# ------------------------------------------------------------
def ovl(output,   # output
        prefix,   # input template
        nf,of,df, # input files n,o,d
        ng,       # number of files in a group
        time=5,   # cluster execution time
        nodes=1): # number of cluster nodes

    # sort groups
    groups=grp(nf,of,df,ng)
    nodes=min(nodes,len(groups))

    # loop over groups
    if(nodes>1): Fork(time=time,ipn=len(groups)/nodes,nodes=nodes)
    for ig in range(len(groups)):
        Plot(output+"-g%04d"%ig,
             [prefix%jg for jg in groups[ig]],
             'Overlay')
        if(nodes>1): Iterate()
    if(nodes>1): Join()

    # join groups
    Plot  (output,[output+"-g%04d"%ig for ig in range(len(groups))],'Overlay')
    Result(output,[output+"-g%04d"%ig for ig in range(len(groups))],'Overlay')

def ovllist(output,   # output
            prefix,   # input template
            list,
            ng,       # number of files in a group
            time=5,   # cluster execution time
            nodes=1): # number of cluster nodes

    # sort groups
    groups=grplist(list,ng)
    nodes=min(nodes,len(groups))

    # loop over groups
    if(nodes>1): Fork(time=time,ipn=len(groups)/nodes,nodes=nodes)
    for ig in range(len(groups)):
        Plot(output+"-g%04d"%ig,
             [prefix%jg for jg in groups[ig]],
             'Overlay')
        if(nodes>1): Iterate()
    if(nodes>1): Join()

    # join groups
    Plot  (output,[output+"-g%04d"%ig for ig in range(len(groups))],'Overlay')
    Result(output,[output+"-g%04d"%ig for ig in range(len(groups))],'Overlay')
    
