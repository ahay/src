from rsf.proj import *
from math import *
import fdmod,pcsutil,wefd

def data(par):

    # ------------------------------------------------------------
    Fetch('vp_marmousi-ii.segy',"marm2")
    Fetch('vs_marmousi-ii.segy',"marm2")
    Fetch('density_marmousi-ii.segy',"marm2")

    # ------------------------------------------------------------
    for file in ('vp','vs','ro'):
        if(file=='ro'):
            ifile='density_marmousi-ii.segy'
            
        else:
            ifile=file+'_marmousi-ii.segy'
            
        Flow(['z'+file,'t'+file,'./s'+file,'./b'+file],ifile,
             '''
             segyread tape=$SOURCE
             tfile=${TARGETS[1]}
             hfile=${TARGETS[2]}
             bfile=${TARGETS[3]}
             ''',stdin=0)
        
        Flow('_'+file,'z'+file,
             '''
             put
             o1=0 d1=0.001249 label1=%(lz)s unit1=%(uz)s
             o2=0 d2=0.001249 label2=%(lx)s unit2=%(ux)s |
             window j1=2 j2=2
             ''' % par)
        
        if(file=='ro'):
            Flow(file+'raw','_'+file,'window n1=%(nz)d n2=%(nx)d min1=%(oz)g min2=%(ox)g | scale rscale=1000000' % par)
        else:
            Flow(file+'raw','_'+file,'window n1=%(nz)d n2=%(nx)d min1=%(oz)g min2=%(ox)g' % par)
                
    # ------------------------------------------------------------
    Flow(  'wmask','vpraw','mask max=1.5 | dd type=float')
#    Result('wmask',fdmod.cgrey('allpos=y',par))
    
    Flow('rx','vpraw','math output="1.0e6+1.5e6*(input-1.5)/3" ')
    Flow('ro','roraw','math output=1')
    
    Flow('vp','vpraw','smooth rect1=35 rect2=35 repeat=5')
    Flow('vs','vp wmask','scale rscale=0.5 | math w=${SOURCES[1]} output="input*(1-w)"')

    # velocity ratio at cig location x
    Flow('vratio1_1','vp vp','add mode=d ${SOURCES[1]}');
    Flow('vratio1_2','vp vs','add mode=d ${SOURCES[1]}');
    Flow('vratio2_1','vs vp','add mode=d ${SOURCES[1]}');
    Flow('vratio2_2','vs vs','add mode=d ${SOURCES[1]}');
    
    Flow('vratio','vratio1_1 vratio1_2 vratio2_1 vratio2_2',
         '''
         cat axis=3 space=n ${SOURCES[0:4]}
         ''',stdin=0)
    
def mask(mask,xsou,tmin,tmax,par):
    
    dipline1(mask+'ml',
                    0.15+tmin,par['xmin'],
                    0.15,xsou,
                    0,1,
                    par['nt'],par['ot'],par['dt'],
                    par['nx'],par['ox'],par['dx'])
    dipline1(mask+'mr',
                    0.15,xsou,
                    0.15+tmax,par['xmax'],
                    0,1,
                    par['nt'],par['ot'],par['dt'],
                    par['nx'],par['ox'],par['dx'])

    Flow(mask,[mask+'ml',mask+'mr'],
         '''
         spike nsp=1 mag=1.0
         n1=%(nx)d o1=%(ox)g d1=%(dx)g k1=%(ltap)d l1=%(rtap)d
         n2=%(nt)d o2=%(ot)g d2=%(dt)g |
         smooth rect1=100 repeat=1 |
         scale axis=123 |
         transp |
         add mode=p ${SOURCES[0]} |
         add mode=p ${SOURCES[1]} |
         transp |
         smooth rect2=100 repeat=3 |
         put label1=x label2=t unit1=km unit2=s |
         spray axis=3 n=2 o=0 d=1 |
         transp plane=23
         ''' % par)
    Result(mask,
           'window n2=1 | transp|' + fdmod.dgrey('',par))
    
def dip(dip,img,par):
    Flow(  dip,img,'dip rect1=40 rect2=40 order=3 liter=100 verb=y ')
    Result(dip,fdmod.cgrey('color=j wantscalebar=n',par))   

def psang(x,img,dip,vpvs,tag,par):
    
    #dip angle at cig location x
    Flow(  dip+'-one',dip,'window n2=1 min2=%g'%x)

    #vpvs ratio at cig location x
    Flow('vratioPP',vpvs,'window n3=1 f3=0 n2=1 min2=%g'%x)
    Flow('vratioPS',vpvs,'window n3=1 f3=1 n2=1 min2=%g'%x)
    Flow('vratioSP',vpvs,'window n3=1 f3=2 n2=1 min2=%g'%x)
    Flow('vratioSS',vpvs,'window n3=1 f3=3 n2=1 min2=%g'%x)
    
    nhx=200
    nhz=0
    nht=0
    
    wefd.elaps('S'+tag,
               img+tag+'_ds',
               img+tag+'_dr',
               nhx,nhz,nht,
               dip+'-one',x,par)



def dipline1(mod,s1,s2,e1,e2,vi,vt,n1,o1,d1,n2,o2,d2):

    min1=o1
    max1=o1+(n1-1)*d1
    min2=o2
    max2=o2+(n2-1)*d2

    ra = (e1-s1)/(e2-s2)
    vels = "%s,%s,%s" %(vi,vt,vt)
    drvs = "%s,%s" %(tan(ra),tan(ra))

    dim1 = 'd1=%g o1=%g n1=%d' % (d2,o2,n2)
    dim2 = 'd2=%g o2=%g n2=%d' % (d1,o1,n1)

    Flow(mod+'lay2',None,
         '''
         spike nsp=4 mag=%g,%g,%g,%g
               n1=4 n2=1 k1=1,2,3,4 |
         put n1=2 n2=2 |
         spline %s fp=%s
         '''%(min2,min1,max2,max1,dim1,drvs))

    Flow(mod+'lay1',None,
         '''
         spike nsp=4 mag=%g,%g,%g,%g
               n1=4 n2=1 k1=1,2,3,4 |
         put n1=2 n2=2 |
         spline %s fp=%s
         '''%(s2,s1,e2,e1,dim1,drvs))
    

    Flow(    mod+'layers',[mod+'lay1',mod+'lay2'],'cat axis=2 ${SOURCES[1:2]}')
    Flow(mod,mod+'layers',
         '''
         unif2 v00=%s n1=%d d1=%g o1=%g
         ''' % (vels,n1,d1,o1) )
