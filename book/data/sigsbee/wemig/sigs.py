from rsf.proj import *
import fdmod

# ------------------------------------------------------------
# model parameters
def param():
    par = {
        'nx':3201,  'ox':10.000,'dx':0.025,  'lx':'x', 'ux':'km',
        'nz':1201,  'oz':0,     'dz':0.025,  'lz':'z', 'uz':'km',
        'nt':1500,  'ot':0,     'dt':0.008,  'lt':'t', 'ut':'s'
        }
    
    par['ft2km']=0.3048
    
    par['ox']=par['ox']*par['ft2km']
    par['dx']=par['dx']*par['ft2km']
    par['oz']=par['oz']*par['ft2km']
    par['dz']=par['dz']*par['ft2km']
    
    return par

# ------------------------------------------------------------
def getdata(data,par):

    datafile = 'sigsbee2a_nfs.sgy'
    Fetch(datafile,'sigsbee')
    
    Flow([data,data+'-t','./'+data+'-h','./'+data+'-b'],
         datafile,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         ''',stdin=0)

# ------------------------------------------------------------
def getmigvel(velo,par):

    migvelfile = 'sigsbee2a_migvel.sgy'
    Fetch(migvelfile,'sigsbee')

    Flow([velo+'-raw',velo+'-t','./'+velo+'-h','./'+velo+'-b'],
         migvelfile,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         ''',stdin=0)
    
    Flow(velo,
         velo+'-raw',
         '''
         scale rscale=0.001 |
         scale rscale=%g |
         put
         o1=%g d1=%g label1=%s unit1=%s
         o2=%g d2=%g label2=%s unit2=%s
         ''' % (par['ft2km'],
                0.0                ,0.0250*par['ft2km'],par['lz'],par['uz'],
                10.025*par['ft2km'],0.0375*par['ft2km'],par['lx'],par['ux']
                ))

# ------------------------------------------------------------
def getstrvel(velo,par):

    strvelfile = 'sigsbee2a_stratigraphy.sgy'
    Fetch(strvelfile,'sigsbee')

    Flow([velo+'-raw',velo+'-t','./'+velo+'-h','./'+velo+'-b'],
         strvelfile,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         ''',stdin=0)

    Flow(velo,
         velo+'-raw',
         '''
         scale rscale=0.001 |
         scale rscale=%g |
         put 
         o1=%g d1=%g label1=%s unit1=%s
         o2=%g d2=%g label2=%s unit2=%s
         ''' % (par['ft2km'],
                0.0                ,0.0250*par['ft2km'],par['lz'],par['uz'],
                10.000*par['ft2km'],0.0250*par['ft2km'],par['lx'],par['ux']
                ))

# ------------------------------------------------------------
def makeshots(shot,data,par):

    Flow(shot+'-ss',data+'-t','dd type=float | headermath output="10925+fldr*150" | window')
    Flow(shot+'-oo',data+'-t','dd type=float | headermath output="offset"         | window')

    # create sraw(t,o,s): o=full offset
    Flow(shot+'-si',shot+'-ss','math output=input/150')
    Flow(shot+'-oi',shot+'-oo','math output=input/75')
    Flow(shot+'-os',[shot+'-oi',shot+'-si'],
         'cat axis=2 space=n ${SOURCES[1]} | transp | dd type=int')
    Flow(shot+'-raw',[data,shot+'-os'],
         '''
         intbin head=${SOURCES[1]} xkey=0 ykey=1 |
         put d2=0.075 d3=0.150 o3=10.95 label1=t label2=o label3=s
         ''')
    Flow(shot,shot+'-raw',
         '''
         mutter half=false t0=1.0 v0=6 |
         put label1=t unit1=s
         o2=%g d2=%g unit2=%s
         o3=%g d3=%g unit3=%s
         ''' % ( 0.000*par['ft2km'],0.075*par['ft2km'],par['ux'],
                 10.95*par['ft2km'],0.150*par['ft2km'],par['ux']
                 ))
        
# ------------------------------------------------------------
def makecmps(cmps,data,par):

    Flow(cmps+'-ss',data+'-t','dd type=float | headermath output="10925+fldr*150" | window')
    Flow(cmps+'-oo',data+'-t','dd type=float | headermath output="offset"         | window')

    Flow(cmps+'-mm',[cmps+'-oo',cmps+'-ss'],
         'math o=${SOURCES[0]} s=${SOURCES[1]} output=s+o/2-10925')
    Flow(cmps+'-hh',[cmps+'-oo'     ],
         'math o=${SOURCES[0]}                 output=o/2')
    
    Flow(cmps+'-mi',cmps+'-mm','math output=input/37.5')
    Flow(cmps+'-hi',cmps+'-hh','math output=input/37.5')
    Flow(cmps+'-mh',[cmps+'-hi',cmps+'-mi'],
         'cat axis=2 space=n ${SOURCES[1]} | transp | dd type=int')
    Flow(cmps+'-raw',[data,cmps+'-mh'],
         '''
         intbin head=${SOURCES[1]} xkey=0 ykey=1 |
         window j2=4 j3=4 n3=500 |
         put d2=0.150 d3=0.150 o3=10.925 label1=t label2=h label3=m
         ''')
    Flow(cmps,cmps+'-raw',
         '''
         mutter half=true t0=0.25 v0=2.250 |
         pad n2out=100 |
         put label1=t unit1=s
         o2=%g d2=%g unit2=%s
         o3=%g d3=%g unit3=%s
         ''' % ( 0.0000*par['ft2km'],0.150*par['ft2km'],par['ux'],
                 10.925*par['ft2km'],0.150*par['ft2km'],par['ux']
                 ))
    

# ------------------------------------------------------------
def symmetrizecmps(symc,cmps,par):
    Flow(symc,cmps,
         '''
         window f2=1 | pad end2=1 |
         reverse which=2 opt=i|
         put o2=%g |
         cat axis=2 space=n ${SOURCES[0]}
         ''' % (-100*0.150*par['ft2km']) )
