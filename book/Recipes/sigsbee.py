try:    from rsf.cluster import *
except: from rsf.proj    import *
import adcig

# ------------------------------------------------------------
# model parameters
def param():
    par = dict(
        nx=3201,  ox=10.000, dx=0.025,  lx='x', ux='km',
        ny=1,     oy=0.000,  dy=0.025,  ly='y', uy='km',
        nz=1201,  oz=0,      dz=0.025,  lz='z', uz='km',
        nt=1500,  ot=0,      dt=0.008,  lt='t', ut='s'
        )

    par['ft2km']=0.3048

    par['ox']=par['ox']*par['ft2km']
    par['dx']=par['dx']*par['ft2km']
    par['oy']=par['oy']*par['ft2km']
    par['dy']=par['dy']*par['ft2km']
    par['oz']=par['oz']*par['ft2km']
    par['dz']=par['dz']*par['ft2km']

    # source coordinates
    par['os']=10.95*par['ft2km']
    par['ds']=0.150*par['ft2km']
    # source index: o=39, d=6, n=500

    par['nzdtm']=244 # redatuming samples through water
    par['nzpad']=143 # padding samples at the bottom

    # all shots parameters
    par['nsall']=500
    par['dsall']=0.04572
    par['osall']=3.33756

    # all receiver parameters
    par['nrall']=348
    par['drall']=0.02286

    # number of nodes
    par['nodes']=8

    # CIG position
    par['xCIG']=9.15

    return par

# ------------------------------------------------------------
def modpar(par):

    par['frq']=10
    par['kt']=120
    par['nt']=12001
    par['dt']=0.001
    par['nb']=150
    par['jsnap']=500
    par['jdata']=1
    par['wweight']=50
    par['wclip']=0.5

# ------------------------------------------------------------
def wempar(par):
    # expand model at the bottom
    par['nz']=par['nz']+par['nzpad']

    par['verb']='y'
    par['eps']=0.1
    par['nrmax']=5
    par['dtmax']=0.00005
    par['tmx']=16

    par['fw']=36
    par['jw']=1
    par['dw']=1/(par['nt']*par['dt'])
    par['kw']=par['nt']/2+1
    par['ow']=par['fw']*par['dw']
    par['nw']=240
    par['eic']='itype=o'

    par['ntpad']=2000

    # define image space
    par['nximg']=1024
    par['jximg']=3
    par['fximg']=65

def fwipar(par):
    par['frq']=10
    par['kt']=120
    par['nt']=12001
    par['dt']=0.001
    par['nb']=150
    par['jsnap']=500
    par['jdata']=1
    par['wweight']=50
    par['wclip']=0.5

    par['nqz']=par['nz']/2
    par['oqz']=par['oz']
    par['dqz']=par['dz']*2

    par['nqx']=par['nx']/2
    par['oqx']=par['ox']
    par['dqx']=par['dx']*2

# ------------------------------------------------------------
def rtmpar(par):
    # expand model at the bottom
    par['nz']=par['nz']+par['nzpad']

    par['frq']=10
    par['kt']=120
    par['nt']=12001
    par['dt']=0.001
    par['nb']=150
    par['jdata']=8
    par['jsnap']=8

    # define image space
    par['nximg']=1024
    par['jximg']=3
    par['fximg']=65

    par['nqz']=par['nz']-par['nzdtm']
    par['oqz']=par['oz']+par['nzdtm']*par['dz']
    par['dqz']=par['dz']

    par['nqx']=par['nximg']
    par['oqx']=par['ox']+par['fximg']*par['dx']
    par['dqx']=par['dx']*par['jximg']

def migpar(par):
    wempar(par)

# ------------------------------------------------------------
def eicpar(par):
    par['nhx']=40
    par['nhy']=0
    par['nhz']=20

    par['nht']=60
    par['dht']=0.008

    adcig.xparam(2*par['nhx']+1,-par['nhx']*par['dx']*par['jximg'],par['dx']*par['jximg'],
                 par['nz']-par['nzdtm'],                 par['oz'],par['dz'],
                 par)
    adcig.tparam(7,
                 2*par['nht']+1,-par['nht']*par['dht'],par['dht'],
                 par['nz']-par['nzdtm'],            par['oz'],par['dz'],
                 par)
    adcig.sparam(5,
                 2*par['nhx']+1,-par['nhx']*par['dx']*par['jximg'],par['dx']*par['jximg'],
                 par['nz'],       par['oz'],                       par['dz'],
                 2*par['nht']+1,-par['nht']*par['dht'],            par['dht'],
                 par)
    adcig.eparam(2,
                 2*par['nhx']+1,-par['nhx']*par['dx']*par['jximg'], par['dx']*par['jximg'],
                 2*par['nhz']+1,-par['nhz']*par['dz'], par['dz'],
                 2*par['nht']+1,-par['nht']*par['dht'],par['dht'],
                 par)

    # slant stack parameters
    par['na']=240
    par['oa']=-60
    par['da']=0.5
    par['ns']=500
    par['os']=-2.5
    par['ds']=0.01

# ------------------------------------------------------------
def hwtpar(par):
    par['ng']=1801
    par['dg']=0.2
    par['og']=-180
    par['lg']='g'
    par['ug']='deg'

# ------------------------------------------------------------
def shotsTWO(par):
    par['fS']=50
    par['jS']=400
    par['nS']=2
    sindex = range(par['fS'],par['fS']+par['nS']*par['jS'],par['jS'])
    par['nodes']=min(par['nodes'],len(sindex))
    return sindex

def shotsWIN(par):
    par['fS']=50
    par['jS']=10
    par['nS']=16
    sindex = range(par['fS'],par['fS']+par['nS']*par['jS'],par['jS'])
    return sindex

def shotsJMP(par):
    par['fS']=0
    par['jS']=4
    par['nS']=120
    sindex = range(par['fS'],par['fS']+par['nS']*par['jS'],par['jS'])
    return sindex

def shotsALL(par):
    par['fS']=0
    par['jS']=1
    par['nS']=480
    sindex = range(par['fS'],par['fS']+par['nS']*par['jS'],par['jS'])
    return sindex

# ------------------------------------------------------------
def getdata(data,par,local=0):

    if(local):
        datafile = 'DATA/sigsbee/sigsbee2a_nfs.sgy'
    else:
        print('local=',local)
        datafile = 'sigsbee2a_nfs.sgy'
        Fetch(datafile,'sigsbee')

    Flow([data,data+'-t','./'+data+'-h','./'+data+'-b'],
        datafile,
        '''
        segyread
        tape=%s
        tfile=${TARGETS[1]}
        hfile=${TARGETS[2]}
        bfile=${TARGETS[3]}
        '''%datafile,stdin=0)

# ------------------------------------------------------------
def getmigvel(velo,par,local=0):

    if(local):
        migvelfile = 'DATA/sigsbee/sigsbee2a_migvel.sgy'
    else:
        migvelfile = 'sigsbee2a_migvel.sgy'
        Fetch(migvelfile,'sigsbee')

        Flow([velo+'-raw',velo+'-t','./'+velo+'-h','./'+velo+'-b'],
            migvelfile,
            '''
            segyread
            tape=%s
            tfile=${TARGETS[1]}
            hfile=${TARGETS[2]}
            bfile=${TARGETS[3]}
            '''%migvelfile,stdin=0)

        Flow(velo,velo+'-raw',
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
def getstrvel(velo,par,local=0):

    if(local):
        strvelfile = 'DATA/sigsbee/sigsbee2a_stratigraphy.sgy'
        Flow([velo+'-raw',velo+'-t','./'+velo+'-h','./'+velo+'-b'],
             strvelfile,
             '''
             segyread
             tape=%s
             tfile=${TARGETS[1]}
             hfile=${TARGETS[2]}
             bfile=${TARGETS[3]}
             '''%strvelfile,stdin=0)
    else:
        strvelfile = 'sigsbee2a_stratigraphy.sgy'
        Fetch(strvelfile,'sigsbee')

        Flow([velo+'-raw',velo+'-t','./'+velo+'-h','./'+velo+'-b'],
        strvelfile,
        	'''
        	segyread
       	 	tape=%s
        	tfile=${TARGETS[1]}
        	hfile=${TARGETS[2]}
        	bfile=${TARGETS[3]}
        	'''%strvelfile,stdin=0)

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
def getreflect(ref,par,local=0):

    if(local):
        reflectfile = 'DATA/sigsbee/sigsbee2a_reflection_coefficients.sgy'
    else:
        reflectfile = 'sigsbee2a_reflection_coefficients.sgy'
        Fetch(reflectfile,'sigsbee')

    Flow([ref+'-raw',ref+'-t','./'+ref+'-h','./'+ref+'-b'],
         None,
         '''
         segyread
         tape=%s
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         '''%reflectfile,stdin=0)

    Flow(ref,
         ref+'-raw',
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
# replace bottom salt and pad with sediment velocity
def replace(vpad,velo,par):

    # padding in z
    par['nreplace']=21
    Flow(velo+'_pad',velo,
         '''
         window n1=1 f1=%d |
         spray axis=1 n=%d |
         smooth rect2=250 repeat=5
         ''' % (par['nz']-par['nreplace'],par['nreplace']+par['nzpad']) )

    Flow(velo+'_sed',velo,'window n1=%d' % (par['nz']-par['nreplace']) )

    Flow(vpad,[velo+'_sed',velo+'_pad'],'cat axis=1 ${SOURCES[1]}')

# ------------------------------------------------------------
# extend bottom salt
def extend(vpad,velo,par):

    # padding in z
    par['nreplace']=21
    Flow(velo+'_cut',velo,
         '''
         window n1=1 f1=%d |
         spray axis=1 n=%d
         ''' % (par['nz']-1,par['nzpad']) )

    Flow(vpad,[velo,velo+'_cut'],'cat axis=1 ${SOURCES[1]}')

# ------------------------------------------------------------
def getrefl(refl,par):

    reffile = 'DATA/sigsbee/sigsbee2a_reflection_coefficients.sgy'

    Flow([refl+'-raw',refl+'-t','./'+refl+'-h','./'+refl+'-b'],
         reffile,
         '''
         segyread
         tape=$SOURCE
         tfile=${TARGETS[1]}
         hfile=${TARGETS[2]}
         bfile=${TARGETS[3]}
         ''',stdin=0)

    Flow(refl,
         refl+'-raw',
         '''
         put
         o1=%g d1=%g label1=%s unit1=%s
         o2=%g d2=%g label2=%s unit2=%s
         ''' % (
                0.0                ,0.0250*par['ft2km'],par['lz'],par['uz'],
                10.000*par['ft2km'],0.0250*par['ft2km'],par['lx'],par['ux']
                ))

# ------------------------------------------------------------
def makeshots(shot,data,par):

    Flow(shot+'-ss',data+'-t',
        'dd type=float | headermath output="10925+fldr*150" | window')
    Flow(shot+'-oo',data+'-t',
        'dd type=float | headermath output="offset"         | window')

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

    Flow(cmps+'-ss',data+'-t',
        'dd type=float | headermath output="10925+fldr*150" | window')
    Flow(cmps+'-oo',data+'-t',
        'dd type=float | headermath output="offset"         | window')

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

# ------------------------------------------------------------
def remap(iout,iinp,imap):
    Flow(iinp+'-transp',iinp,'transp')
    Flow(imap+'-transp',imap,'transp')
    Flow(iout+'-temp',[iinp+'-transp',imap+'-transp'],
         'remap1 pattern=${SOURCES[1]} | transp')
    Flow(iout,[iout+'-temp',imap],
         'remap1 pattern=${SOURCES[1]}')

# ------------------------------------------------------------
def makemask(velo,smask,wmask,lmask,par):
    # salt mask
    Flow(  smask,velo,'mask min=4.499 | dd type=float')

    # water mask
    Flow(  wmask,velo,'mask max=1.5 | dd type=float')

    # sediment mask
    Flow(lmask,[smask,wmask],'add ${SOURCES[1]} | math output="1-input"')

# ------------------------------------------------------------
# low velocity (in sediments only)
def makevlow(vlow,velo,smask,wmask,lmask,coef):
    Flow(vlow,[velo,smask,wmask,lmask],
         '''
         math
         v=${SOURCES[0]}
         s=${SOURCES[1]}
         w=${SOURCES[2]}
         l=${SOURCES[3]}
         output="(%f)*v*l+1.5*w+4.5*s"
         ''' % coef)

# high velocity (in sediments only)
def makevhig(vhig,velo,smask,wmask,lmask,coef):
    Flow(vhig,[velo,smask,wmask,lmask],
         '''
         math
         v=${SOURCES[0]}
         s=${SOURCES[1]}
         w=${SOURCES[2]}
         l=${SOURCES[3]}
         output="(%f)*v*l+1.5*w+4.5*s" |
         add add=-11.5 |
         clip clip=10 |
         add add=11.5
         ''' % coef)

# ------------------------------------------------------------
# density by Gardner's relation
def makedensity(velo,dens,smask,wmask,lmask,par):

    Flow(dens,[velo,smask,wmask,lmask],
         '''
         math
         v=${SOURCES[0]}
         s=${SOURCES[1]}
         w=${SOURCES[2]}
         l=${SOURCES[3]}
         output="0.23*((v*l)*3280.0)^0.25+1.0*w+1.23*s"
         ''')
