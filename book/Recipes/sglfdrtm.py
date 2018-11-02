from rsf.proj import * 
import sys  

# --------------------------------------------------------------------------
# Example of parameters for sglfdrtm
rtmpar = {
    'lnx':  5395,
    'lnz':  956,
    'nx' :  3201,
    'nz' :  956,
    'nt' :  7500,
    'dx' :  12.5,
    'dz' :  12.5,
    'dt' :  0.0016,
    'labelx': "Distance",
    'labelz': "Depth",
    'unitx' : "Km",
    'unitz' : "Km",
    'shtbgn':  2001,
    'shtend':  3595,
    'sintv' : 4,
    'spz'   : 10,
    'gp'    : 0,
    'size'   : 8,
    'frqcut' : 1.0,
    'pml'    : 30,
    #source
    'srcbgn'  : 100,
    'frq'     : 20

    }

# --------------------------------------------------------------------------

def igrey(custom):
    return '''
    grey color=j labelsz=5 titlesz=6 %s
    ''' %(custom)

def grey(custom):
    return '''
    grey labelsz=5 titlesz=6 %s
    ''' %(custom)

def getpar(ipar, idn):
    opar = ipar
    opar['shtnum']   = int((ipar['shtend']-ipar['shtbgn'])/ipar['sintv'])+1
    opar['curid']    = idn
    opar['cursht']   = ipar['shtbgn']  + ipar['sintv']*idn
    opar['curxbgn']  = opar['cursht']  - int(ipar['nx']*0.5)
    opar['curxend']  = opar['curxbgn'] + ipar['nx'] - 1
    opar['spx']      = int(ipar['nx']*0.5)
    opar['bd']       = int(opar['size']*0.5+0.5) + opar['pml']
    checkpar(opar)
    return opar


def checkpar(par):
    if par['cursht'] > par['shtend'] or par['cursht'] < par['shtbgn'] :
        print(" >>>> current source error cursht=%d " %(par['cursht'])) 
        sys.exit()
    
    if par['curxbgn'] < 0 or par['curxend'] > par['lnx'] :
        print(" >>>> curxbgn error: curxbgn=%d, curxend=%d" %(par['curxbgn'],  par['curxend'])) 
        sys.exit()
    
def printpar(par):
    keylist = list(par.keys())
    keylist.sort()
    print("{")
    for key in keylist:
        print('  {0:<10}'.format(key)+":  %s" %par[key]);
    print("}")
           
# --------------------------------------------------------------------------     
def splitmodel(Fout, Fin, par):
    Flow(Fout,Fin,
         '''
         window n1=%(nz)d n2=%(nx)d f2=%(curxbgn)d |
         put label1=%(labelz)s unit1=%(unitz)s
             label2=%(labelx)s unit2=%(unitz)s 
         '''%par)
    
def sglfdrtm(Fimg1, Fimg2, Fsrc, Ffvel, Ffden, Fbvel, Fbden, par, surfix):
    
    Gx = 'Gx%s'    %surfix
    Gz = 'Gz%s'    %surfix
    sxx ='sxx%s'   %surfix 
    sxz ='sxz%s'   %surfix
    szx = 'szx%s'  %surfix
    szz = 'szz%s'  %surfix
    Frcd = 'record%s' %surfix 
    
    Ftmpwf =  'tmpwf%s'  %surfix
    Ftmpbwf = 'tmpbwf%s' %surfix

#  -------------------------------------------------------------------
#  comment this part when use pscons 
#  
    for m in [Ffden, Ffvel,Fbden, Fbvel]:
        pml  = m+'_pml'
        Flow(pml,m,
             '''
             expand left=%(bd)d right=%(bd)d 
                    top=%(bd)d  bottom=%(bd)d
             '''%par)
#  -------------------------------------------------------------------
    Ffdenpml = Ffden+'_pml'
    Ffvelpml = Ffvel+'_pml'
    Fbdenpml = Fbden+'_pml'
    Fbvelpml = Fbvel+'_pml'
    
#  -------------------------------------------------------------------
#  comment this part when use pscons 
#
    Flow([Gx,sxx,sxz],Ffvelpml,
         '''
         sfsglfdcx2_7 dt=%(dt)g eps=0.00001 npk=50 
                      size=%(size)d sx=${TARGETS[1]} sz=${TARGETS[2]}
                      wavnumcut=%(frqcut)g
         '''%par)
#  -------------------------------------------------------------------
#  comment this part when use pscons 
#
    Flow([Gz,szx,szz],Ffvelpml,
         '''
         sfsglfdcz2_7 dt=%(dt)g eps=0.00001 npk=50 
                      size=%(size)d sx=${TARGETS[1]} sz=${TARGETS[2]}
                      wavnumcut=%(frqcut)g
         '''%par)
#  -------------------------------------------------------------------    
    Flow([Fimg1, Fimg2, Frcd],[Fsrc,Ffvelpml,Ffdenpml,Gx,sxx,sxz,Gz,szx,szz,Fbvelpml,Fbdenpml],
       '''sfsglfdrtm2 img2=${TARGETS[1]} rec=${TARGETS[2]} 
          fvel=${SOURCES[1]} fden=${SOURCES[2]}
          bvel=${SOURCES[9]} bden=${SOURCES[10]}
          Gx=${SOURCES[3]} sxx=${SOURCES[4]} sxz=${SOURCES[5]}
          Gz=${SOURCES[6]} szx=${SOURCES[7]} szz=${SOURCES[8]}
          freesurface=n  verb=y decay=n 
          spx=%(spx)g spz=%(spz)g pmlsize=%(pml)d snapinter=10 
          srcdecay=y  gp=%(gp)g srctrunc=0.2 pmld0=200
          wantrecord=n
       '''%par)

#  -------------------------------------------------------------------
#  test for pscons
def sglfdrtm_test(Fimg1, Fimg2, Fsrc, Fvel, Fden, par, surfix):
    
    Gx   = 'Gx%s'    %surfix
    Gz   = 'Gz%s'    %surfix
    sxx  = 'sxx%s'   %surfix 
    sxz  = 'sxz%s'   %surfix
    szx  = 'szx%s'   %surfix
    szz  = 'szz%s'   %surfix
    Frcd = 'record%s'%surfix 
    
    Ftmpwf  =  'tmpwf%s'  %surfix
    Ftmpbwf =  'tmpbwf%s' %surfix

    Fdenpml = Fden+'_pml'
    Fvelpml = Fvel+'_pml'

    Flow([Fimg1, Fimg2, Frcd, Fvelpml,Fdenpml,Gx,sxx,sxz,Gz,szx,szz],[Fsrc,Fvel,Fden],
         '''
         sfexpand   <${SOURCES[1]}
                    left=%(bd)d right=%(bd)d 
                    top=%(bd)d  bottom=%(bd)d
                    >${TARGETS[3]}  &&
         
         sfexpand   <${SOURCES[2]}
                    left=%(bd)d right=%(bd)d 
                    top=%(bd)d  bottom=%(bd)d
                    >${TARGETS[4]}  &&
         
         sfsglfdcx2_7  <${TARGETS[3]} 
                       dt=%(dt)g eps=0.00001 npk=50 
                       size=%(size)d sx=${TARGETS[6]} sz=${TARGETS[7]}
                       wavnumcut=%(frqcut)g
                       >${TARGETS[5]}  &&

         sfsglfdcz2_7  <${TARGETS[3]}
                       dt=%(dt)g eps=0.00001 npk=50 
                       size=%(size)d sx=${TARGETS[9]} sz=${TARGETS[10]}
                       wavnumcut=%(frqcut)g
                       >${TARGETS[8]}  &&

         sfsglfdrtm2   <${SOURCES[0]}
                       img2=${TARGETS[1]} rec=${TARGETS[2]} 
                       fvel=${TARGETS[3]} fden=${TARGETS[4]}
                       bvel=${TARGETS[3]} bden=${TARGETS[3]}
                       Gx=${TARGETS[5]} sxx=${TARGETS[6]} sxz=${TARGETS[7]}
                       Gz=${TARGETS[8]} szx=${TARGETS[9]} szz=${TARGETS[10]}
                       freesurface=n  verb=y decay=n 
                       spx=%(spx)g spz=%(spz)g pmlsize=%(pml)d snapinter=10 
                       srcdecay=y  gp=%(gp)g srctrunc=0.2 pmld0=200
                       >${TARGETS[0]}
         
         '''%par,stdin=0, stdout=0)

    

        



