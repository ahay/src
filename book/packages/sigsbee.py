from rsfproj import *
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

    par['nb']=250

    # source coordinates
    par['os']=10.95*par['ft2km']
    par['ds']=0.150*par['ft2km']
    # source index: o=39, d=6, n=500
    
    # receiver coordinates
    par['or']=10.95*par['ft2km']
    par['dr']=0.075*par['ft2km']

    par['nzdtm']=244 # number of redatuming steps through water
    par['nzpad']=143

    # all shots parameters
    par['nsall']=500
    par['dsall']=0.04572
    par['osall']=3.33756

    return par

# ------------------------------------------------------------
def modpar(par):

    par['kt']=100
    par['nt']=12001
    par['dt']=0.001
    par['nb']=150
    par['jsnap']=1000
    par['jdata']=1
    par['wweight']=50
    par['wclip']=0.5
 
# ------------------------------------------------------------   
def migpar(par):

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
    
# ------------------------------------------------------------
def getdata(data,par):

    datafile = 'data/sigsbee/sigsbee2a_nfs.sgy'
    #Fetch(data,'sigsbee')
    
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

    migvelfile = 'data/sigsbee/sigsbee2a_migvel.sgy'
    #Fetch(velo,'sigsbee')

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

    strvelfile = 'data/sigsbee/sigsbee2a_stratigraphy.sgy'
    #Fetch(velo,'sigsbee')

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
         o2=%g d2=%g label2=%s unit2=%s |
	 window n1=%d n2=%d
         ''' % (par['ft2km'],
                0.0                ,0.0250*par['ft2km'],par['lz'],par['uz'],
                10.000*par['ft2km'],0.0250*par['ft2km'],par['lx'],par['ux'],
		par['nz'],par['nx']
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
def getrefl(refl,par):

    reffile = 'data/sigsbee/sigsbee2a_reflection_coefficients.sgy'

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
    Result(smask,fdmod.cgrey('allpos=y',par))
    
    # water mask
    Flow(  wmask,velo,'mask max=1.5 | dd type=float')
    Result(wmask,fdmod.cgrey('allpos=y',par))

    # sediment mask
    Flow(lmask,[smask,wmask],'add ${SOURCES[1]} | math output="1-input"')
    Result(lmask,fdmod.cgrey('allpos=y',par))

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
    
