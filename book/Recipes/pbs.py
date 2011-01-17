try:
    from rsf.cluster import *
except:
    from rsf.proj import *

''' 
A series of functions that are used by the sfpbswrap program for
individual clusters.
'''

def setSConsign(env,suffix):
    try:     
        import dbhash
        env.SConsignFile(env.path+ '.sconsign-%s.dbhash' % suffix,dbhash)
    except:      
        try:     
            import gdbm
            env.SConsignFile(env.path+ '.sconsign-%s.gdbm' % suffix ,gdbm)
        except: 
            env.SConsignFile(env.path+ '.sconsign-%s' % suffix)

