from rsf.proj import *
import sigs,fdmod,pot


def inimodel(ss,rr,par):

    
    Result('vwin',fdmod.cgrey('color=j allpos=y wantscalebar=y',par))
    
    # ------------------------------------------------------------
    sigs.saltmask('mask','vwin',par)
    Result('mask',fdmod.cgrey('wantscalebar=y',par))
    
    Flow('vo',['vwin','mask'],
         'math m=${SOURCES[1]} output="input*(1-m)+m*0.85*14.76*%g"' % par['ft2km'])

    # ------------------------------------------------------------
    Flow('vp','vwin','scale rscale=1')
    Flow('nu','vp','math output=0.0')
    
    # vp/vs
    Flow('vsvpratio','vwin',
         '''
         math output="0.40+0.25*(x1-2.)/9."
         ''')
    Flow('vs','vp vsvpratio mask',
         '''
         math r=${SOURCES[1]} output="input*r" |
         math m=${SOURCES[2]} output="input*(1-m)+m*8*%g"
         ''' % par['ft2km'])

    # epsilon and delta
    Flow(  'amask','vp','mask min=2.0 max=4 | dd type=float')
    Result('amask',fdmod.cgrey('allpos=y',par))
    Flow(  'bmask','vp','mask min=2.5 max=3.0 | dd type=float')
    Result('bmask',fdmod.cgrey('allpos=y',par))    

    Flow(  'cmask','vp','mask min=2.5 max=4 | dd type=float')
    Result('cmask',fdmod.cgrey('allpos=y',par))
    
    Flow(  'denmask','vp','mask max=4 | dd type=float')
    Result('denmask',fdmod.cgrey('allpos=y',par))


#    Flow('ro','vp','math output="-input*0.27+3.3" ')
    # ro
    Flow('ro',['mask','vp'],
         '''
         math output="-a "      a=${SOURCES[0]}  |
         math output="input+v*0.8-.2"   v=${SOURCES[1]}
         ''')
    # epsilon
    Flow('epsilon',['amask','bmask','vp'],
         '''
         math output="a*0.28+b*0.2"
         a=${SOURCES[0]} b=${SOURCES[1]} |
         math output="input*v/3.0"
         v=${SOURCES[2]}
         ''')

    # delta
    Flow('delta',['amask','bmask','cmask','vp'],
         '''
         math output="a*0.0+b*0.05+c*0.05"
         a=${SOURCES[0]} b=${SOURCES[1]} c=${SOURCES[2]} |
         math output="input*v/3.0"
         v=${SOURCES[3]}
         ''')
    labelattr='labelsz=11 wantscalebar=y bartype=h yll=1.5 screenht=9.5 wherebarlabel=top wherebartics=top wherebar=top '
    Plot('vp',
         fdmod.cgrey('allpos=y bias=1.43 color=j barunit="km/s" barlabel="V\_P\^" '
                     +labelattr,par))
    Plot('vs',
         fdmod.cgrey('allpos=y bias=0.72 color=j barunit="km/s" barlabel="V\_S\^" '
                     +labelattr,par))
    Plot('ro',
         fdmod.cgrey('allpos=y bias=1.0 color=j barunit="\F4 g/cm\^3" barlabel="\F9 r" '
                     +labelattr,par))
    Plot('epsilon',
         fdmod.cgrey('allpos=y color=iC color=j formatbar=%4.2f  barlabel="\s150 \F9 e" '
                     +labelattr,par))
    Plot('delta',  
         fdmod.cgrey('allpos=y color=iC color=j formatbar=%4.2f barlabel="\s150 \F9 d" '
                     +labelattr,par))



    Plot('vpvsratio','vsvpratio','math output="1/input" |'+fdmod.cgrey('allpos=y color=iC color=j',par))
    for i in (['vp','vs','ro','epsilon','delta','vpvsratio']):
        Result(i,[i],'Overlay')
    
def wfom(wom,wfldp,wflds,velo,vmean,name1,name2,axis,custom,par):
    
#    if(not par.has_key('wweight')): par['wweight']=10
#    if(not par.has_key('wclip')):   par['wclip']=1.0

    print(par['wweight'])

    for mode in [wfldp,wflds]:
        Flow(mode+'t',mode,'scale axis=123 |scale rscale=5')
        Flow(mode+'t2',[velo,mode+'t'],
             '''
             add add=-%g |
             scale axis=123 |
             math w=${SOURCES[1]} output="input+%g*w"
             ''' %(vmean, par['wweight']) )

    pot.cliptogether(wom,wfldp+'t2',wflds+'t2',name1,name2,axis,custom,par)

def wfom3(wom,wfldp,wfldsv,wfldsh,velo,vmean,name1,name2,name3,axis,custom,par):

    for mode in [wfldp,wfldsv,wfldsh]:
        Flow(mode+'t',mode,'scale axis=123 |scale rscale=5')
        Flow(mode+'t2',[velo,mode+'t'],
             '''
             add add=-%g |
             scale axis=123 |
             math w=${SOURCES[1]} output="input+%g*w"
             ''' %(vmean, par['wweight']) )
        
    pot.cliptogether3(wom,wfldp+'t2',wfldsv+'t2',wfldsh+'t2',name1,name2,name3,axis,custom,par)    





def weight(nref,par,param): 

    Flow('re','R epsilon','add mode=p ${SOURCES[1]}')

    for i in range(0,nref):        
        ref="%01d"%i

                
#        Flow('eps'+ref,'epsilon',
#             '''
#             window n1=1 n2=1 min1=%f min2=%f |
#             spray n=%d axis=1 | spray n=%d  axis=2
#             '''%(param['z'+ref],param['x'+ref],par['nz'],par['nx']) 
#             )
#        Flow('R'+ref,'R',
#             '''
#             window n1=1 n2=1 min1=%f min2=%f |
#             spray n=%d axis=1 | spray n=%d  axis=2
#             '''%(param['z'+ref],param['x'+ref],par['nz'],par['nx']) 
#             )
        Flow('re'+ref,'re',
             '''
             window n1=1 n2=1 min1=%f min2=%f |
             spray n=%d axis=1 | spray n=%d  axis=2
             '''%(param['z'+ref],param['x'+ref],par['nz'],par['nx']) 
             )
#        Flow('w'+ref,['epsilon','R','eps'+ref,'R'+ref],
#             '''
#             math output="sqrt((a-a0)*(a-a0)+(b-b0)*(b-b0))"
#             a=${SOURCES[0]} a0=${SOURCES[1]}
#             b=${SOURCES[2]} b0=${SOURCES[3]}
#             ''')
        Flow('w'+ref,['re','re'+ref],
             '''
             math output="abs(a-a0)"
             a=${SOURCES[0]} a0=${SOURCES[1]}
             ''')
