try:
    from rsf.cluster import *
except:
    from rsf.proj import *
import commands,os


'''
    Module designed to get interaction with useful 
    shell m8r utilities as: sfattr, sfdisfil, sfin, etc
    
    It should be extended to take advantage of other unix
    utilities (ls, tail, grep, etc)

    Esteban Diaz,
    Center for Wave Phenomena
'''


def sfaxa(fin,dim=1): 
    '''
    return a vector with the geometry along dimension dim
    axa=['n1','d1','o1'] 
    '''
    var=commands.getoutput('sfin %s.rsf'%fin+'|grep n%d | sed -e \"s/=/ /g\"'%dim)
    tmp=var.split()
    axa=[float(tmp[1]),float(tmp[3]),float(tmp[5])] 
    return axa

def sfdisfil(fin):
    '''
    return the printout of sfdisfil
    '''
    dis=commands.getoutput('sfdisfil< %s.rsf'%fin)
    return dis

def sfattr(fin,attr='min'):
    '''
    return the a user defined attribute from file fin,
    posibles attr values are:
    rms, mean, 2-norm, variance, std dev, max, min     
    '''
    tmp=commands.getoutput('sfattr< %s.rsf'%fin+'|grep %s'%attr)
    tmp2=tmp.split()
    attr=tmp2[2]
    return attr


def rmbin(fin):
    '''
    remove the binary from a rsf file
    (which scons does not have a clue about it)
    '''
    var=commands.getoutput('sfin %s.rsf'%fin+'|grep in | sed -e "s/\\"/ /g"')
    tmp=var.split()
    cmd='rm -f %s'%tmp[2]
    os.system(cmd)

def rmtarget(fin):
    '''
    remove fin file
    '''
    cmd='rm -f %s'%fin
    os.system(cmd)

