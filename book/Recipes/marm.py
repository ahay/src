from rsf.proj import *
import fdmod,wefd,pcsutil,zomig,spmig,pplot,fdd,stiff

def param():
    
    par = {
        'nt':2001,'ot':0, 'dt':0.0005, 'lt':'t', 'ut':'s',
        'nx':13601, 'ox':8.00, 'dx':0.001249, 'lx':'x', 'ux':'km',
        'ny':100, 'oy':8.00, 'dy':0.001249, 'ly':'y', 'uy':'km',
        'nz':2801, 'oz':0.50, 'dz':0.001249, 'lz':'z', 'uz':'km',
        'kt':200,
        'frq':15,
        'jsnap':250,
        'jdata':8,
        'nb':20,
        'ng':400,'dg':0.2,'og':-40,
        'wclip':0.5,
        'wweight':5,
        'iweight':25
        }
    par['jz']=4
    par['jx']=4
    par['jy']=4
   
    par['xmax']=14.0
    par['zmax']=3.5
    par['ymax']=9.0
    
    par['dz']=par['dz'] * par['jz']
    par['dx']=par['dx'] * par['jx']
    par['dy']=par['dy'] * par['jy']
    par['nz']=int( (par['zmax'] -par['oz'])/par['dz'] )+1
    par['nx']=int( (par['xmax'] -par['ox'])/par['dx'] )+1
    par['ny']=int( (par['ymax'] -par['oy'])/par['dy'] )+2
   


    fdmod.param(par)
    par['labelattr']=par['labelattr']+'''
        titlesz=12 labelsz=6 bartype=h  font=4 barlabelsz=6 barlabelfat=3
        wantaxis=y wantlabel=y wanttitle=n  parallel2=n  d2num=1 o2num=1 n2tic=3
        wantaxis=y wantlabel=y wanttitle=n 
        '''
    par['xsou']=9.5
    par['zsou']=par['oz']
    par['zrec']=par['oz']
    
    par['ftap']=100
    par['ntap']=par['nx']-par['ftap']

    par['kbem']=680
    par['lbem']=721

    par['nqz']=par['nz']
    par['nqx']=par['nx']
    par['dqz']=par['dz']
    par['dqx']=par['dx']
    par['oqz']=par['oz']
    par['oqx']=par['ox']

    par['fw']=par['nt']*par['dt']
    par['kw']=par['nt']/2+1
    par['dw']=1/(par['nt']*par['dt'])
    par['nw']=150
    par['jw']=1
    par['ow']=5*par['dw']
    par['verb']='y'
    par['nrmax']=3
    par['dtmax']=0.0001

    par['models']=['0','1']

    par['jtframe'] = 50
    par['ntframe'] = par['nt']/par['jdata']/par['jtframe']

    par['jzframe'] = 100
    par['nzframe'] = par['nz']/par['jzframe']

    par['jwframe'] = 10
    par['nwframe'] = par['nw']/par['jwframe']
    
  
    return par

# ------------------------------------------------------------
def grabdata(par):
    
    # ------------------------------------------------------------
    Fetch('vp_marmousi-ii.segy',"marm2")
    Fetch('density_marmousi-ii.segy',"marm2")
    
    # ------------------------------------------------------------
    for file in ('vp','rx','vs'):
        if(file=='vp'):
            ifile='vp_marmousi-ii.segy'
        elif(file=='vs'):
            ifile='vs_marmousi-ii.segy'
        else:
            ifile='density_marmousi-ii.segy'
            
        Flow(['z'+file,'t'+file,'./s'+file,'./b'+file],ifile,
             '''
             segyread tape=$SOURCE
             tfile=${TARGETS[1]}
             hfile=${TARGETS[2]}
             bfile=${TARGETS[3]}
             ''',stdin=0)
     
        Flow(file+'win','z'+file,
             '''
             put
             o1=0 d1=0.001249
             o2=0 d2=0.001249 |
             window  min1=%(oz)g  j1=%(jz)d
                     min2=%(ox)g  j2=%(jx)d n2=%(nx)d 
             ''' % par)
       
# ------------------------------------------------------------
def model(vp,vs,rx,epsilon,delta,nu,par):
    label=''
    barattr=''' xll=2 yll=1.3 
         wantscalebar=y bartype=h
         wherebartic=top wherebarlabel=top barlabelsz=8
         labelsz=8 color=j'''
    Result(vp,fdmod.cgrey('bias=1.6 min1=.5  barlabel="V\_P0\^ (km/s)" '+barattr,par))
    Result(vs,fdmod.cgrey('bias=0.0 min1=.5  barlabel="V\_S0\^ (km/s)" '+barattr,par))
    Result(rx,fdmod.cgrey('bias=1.72 min1=.5 barlabel="\s140 \F9 r \F4 (g/cm\^3) " '+barattr,par))
    Result(epsilon,fdmod.cgrey('bias=.1 min1=.5 barlabel="\s150 \F9 e " formatbar=%4.2f'+barattr,par))
    Result(delta,fdmod.cgrey('bias=0.1 min1=.5  barlabel="\s150 \F9 d " formatbar=%4.2f'+barattr,par))
    Result(nu,fdmod.cgrey('bias=0 color=e allpos=n min1=.5 formatbar=%3.0f barlabel="\s150 \F10 n \s100 (\^o\_)" '+barattr+' color=e',par))
    
## ------------------------------------------------------------ 

   
def dip(dip,img,par):

   Flow(dip+'-tmp',img,'dip rect1=5 rect2=5 order=3 liter=100 verb=y ')

   Flow(dip+'-top',dip+'-tmp',
        '''window f1=0 n1=450 |
           add add=10 | clip clip=10 | add add=-10
        ''')
   Flow(dip+'-bot',dip+'-tmp',
        '''
        window f1=450
        ''')
   Flow(dip,[dip+'-top',dip+'-bot'],'''cat axis=1 space=n ${SOURCES[1]} |
             math output=" atan(input)*180/3.14159265" ''')
#   Result(dip,fdmod.cgrey('color=E allpos=n',par))





def xplot(epsilon,delta,nu):

    # ------------------------------------------------------------
    minn=-60; maxn=60;
    mine=0.10; maxe=0.35;
    mind=-0.10;maxd=0.25;

    xplotformat=''' xll=3.5 yll=2.5 title= symbol=.
    parallel2=n min1=%g max1=%g min2=%g max2=%g  
    labelsz=12 screenratio=1'''

    Plot('enplot',[epsilon,nu],
         '''cmplx ${SOURCES[1]}  |   
            put n1=360000 n2=1 n3=1 | 
            window j1=20 |
            graph label1="\s150 \F9 e" unit1=  
                    label2="\F9 \s150 n \F4 \s100  (\^o\_)" '''+xplotformat
         %(mine,maxe,minn,maxn) )
    Plot('edplot',[epsilon, delta],
         '''cmplx ${SOURCES[1]}  |   
            put n1=360000 n2=1 n3=1 | 
            window j1=20 |
            graph label1="\F9 \s150 e" unit1= 
                    label2="\F9 \s150 d" '''+xplotformat
         %(mine,maxe,mind,maxd) )

    Plot('ndplot',[nu,delta],
         '''cmplx ${SOURCES[1]}  |   
            put n1=360000 n2=1 n3=1 | 
            window j1=20 |
            graph label1="\F9 \s150 n \F4 \s100  (\^o\_)"  unit1= 
                    label2="\F9  \s150 d" l'''+xplotformat
         %(minn,maxn,mind,maxd) )
    pplot.p2x2('xplot','enplot','','edplot','ndplot',0.5,0.5,-10,-12)

def xref(ref):

    minn=-60; maxn=60;
    mine=0.10; maxe=0.35;
    mind=-0.10;maxd=0.25;
    xrefformat=''' xll=3.5 yll=2.5  title= symbol=x symbolsz=10 title= 
                   parallel2=n min1=%g max1=%g min2=%g max2=%g 
                   labelsz=12 screenratio=1'''
    Flow('epref',ref,'window  n2=1 f2=0')
    Flow('deref',ref,'window  n2=1 f2=1')
    Flow('nuref',ref,'window  n2=1 f2=2')
    Plot('enref','epref nuref',
       '''cmplx ${SOURCES[1]}  | 
          graph label1="\s150 \F9 e" unit1=  
                label2="\F9 \s150 n \F4 \s100  (\^o\_)" '''+xrefformat
     %(mine,maxe,minn,maxn) )
    Plot('edref','epref deref',
       '''cmplx ${SOURCES[1]} | 
          graph label1="\F9 \s150 e" unit1= 
                label2="\F9 \s150 d" '''+xrefformat
     %(mine,maxe,mind,maxd) )
    Plot('ndref','nuref deref',
       '''cmplx ${SOURCES[1]}| 
          graph label1="\F9 \s150 n \F4 \s100  (\^o\_)"  unit1= 
                label2="\F9  \s150 d" l'''+xrefformat
     %(minn,maxn,mind,maxd) )
    pplot.p2x2('xref','enref','','edref','ndref',0.5,0.5,-10,-12)



def density(ref,den,x,y,z,w,custom1,custom2):
    Flow(den,[x, y, z, w],
     ''' density n1=30 n2=30 n3=5 n4=1
         inY=${SOURCES[1]} inZ=${SOURCES[2]} inW=${SOURCES[3]} 
         verb=y '''+custom1)
    
    Flow(ref,den,'sortdensity '+custom2)


# ------------------------------------------------------------
def AnalyInterp2d(suffix,uAz,uAx,cAref,nepdel,nnu,epsilon,delta,nu,par):   
    
    listPnu=[]
    listSnu=[]
    for j in range(0,nnu):
        refnu="%01d"%j
        nux=nu[j]
        print nux
        listP=[]
        listS=[]
        Flow('nu'+refnu,'','spike n1=1 mag=%f'%nux)
        for i in range(0,nepdel): 
            ref="%01d"%i        
            ep = epsilon[i]
            de = delta[i]
            ca = stiff.tti2d_point(3, 1.5, 2.4, ep, de, nux)
            Flow('cA'+ref+'_nu'+refnu,'',
                 'spike n3=6 nsp=6 k3=1,2,3,4,5,6 mag=%f,%f,%f,%f,%f,%f n1=1 n2=1 '
                 %(ca[0],ca[1],ca[2],ca[3],ca[4],ca[5]) )

            fdd.separatorD('dzK'+ref+'_nu'+refnu,'dxK'+ref+'_nu'+refnu,
                          'spkk',
                          'cA'+ref+'_nu'+refnu,
                          'y','k','gauss',0.866 ,8, 25, 25,par)

            listP.append('PP'+ref+'_nu'+refnu)
            listS.append('SS'+ref+'_nu'+refnu)
            fdd.SepK('P'+ref+'_nu'+refnu,'S'+ref+'_nu'+refnu,
                     'uAz','uAx',
                     'dzK'+ref+'_nu'+refnu+'-tmp','dxK'+ref+'_nu'+refnu+'-tmp',
                     par)
            Flow('PP'+ref+'_nu'+refnu,['wtepdel','P'+ref+'_nu'+refnu],'window n3=1 f3=%d |add ${SOURCES[1]} mode=p'%(i) )
            Flow('SS'+ref+'_nu'+refnu,['wtepdel','S'+ref+'_nu'+refnu],'window n3=1 f3=%d |add ${SOURCES[1]} mode=p'%(i) )
        Flow('PP_nu'+refnu,listP,'add ${SOURCES[1:%d]}'%nepdel)    
        Flow('SS_nu'+refnu,listS,'add ${SOURCES[1:%d]}'%nepdel)    


        listPnu.append('sepp_nu'+refnu+'_P_wtd')
        listSnu.append('sepp_nu'+refnu+'_S_wtd')

        Flow('sepp_nu'+refnu+'_P_wtd',['wtnu','PP_nu'+refnu],
             'window n3=1 f3=%d |add ${SOURCES[1]} mode=p'%(j) )
        Flow('sepp_nu'+refnu+'_S_wtd',['wtnu','SS_nu'+refnu],
             'window n3=1 f3=%d |add ${SOURCES[1]} mode=p'%(j) )

    Flow('PP-'+suffix,listPnu,'add ${SOURCES[1:%d]}'%nepdel)    
    Flow('SS-'+suffix,listSnu,'add ${SOURCES[1:%d]}'%nepdel)    

    
# ------------------------------------------------------------

def interpSep2d(wtlst,uAz,uAx,cAref,n,par):   

    for j in range(n):

         ref="-%02d"%j     
         Flow('cA'+ref,'cAref','window n1=1 f1=%d squeeze=n'%j)
         Flow('nu'+ref,'nuref','window n1=1 f1=%d'%j)
         fdd.separatorD('dzK'+ref,'dxK'+ref,
                        'spkk',
                        'cA'+ref,
                        'y','k','gauss', 0.866, 8, 25, 25, par)
         fdd.SepK('P'+ref,'S'+ref,
                  uAz,uAx,
                  'dzK'+ref+'-tmp','dxK'+ref+'-tmp',
                  par)

         
         for wt in wtlst:
             Flow('PP'+ref+'-'+wt,['wt'+wt,'P'+ref],
                  'window n3=1 f3=%d | add ${SOURCES[1]} mode=p'%j )
             Flow('SS'+ref+'-'+wt,['wt'+wt,'S'+ref],
                  'window n3=1 f3=%d | add ${SOURCES[1]} mode=p'%j )


    for wt in wtlst:
        Flow('PP-'+wt,
             map(lambda x: 'PP-%02d-' % x +wt, range(n)),
             'add ${SOURCES[1:%d]}'%n )
        Flow('SS-'+wt,
             map(lambda x: 'SS-%02d-' % x +wt, range(n)),
                 'add ${SOURCES[1:%d]}'%n )



def interpSep3d(wt,uAz,uAx,uAy,cAref,n,par,nuu,alpha):   
    print wt
    for j in range(n):

         ref=wt+"-%02d"%j     
         Flow('cA'+ref,'cAref','window n1=1 f1=%d squeeze=y'%j)
         Flow('nu'+ref,'nuref','window n1=1 f1=%d'%j)
         Flow(['dzK'+ref,'dxK'+ref,'dyK'+ref],
	     ['spkk3', 'cA'+ref,'Code/EDERIV3TTI_v.x'],
	     '''
             ${SOURCES[2]} verb=y stat=y domain=k ompnth=8
             ccc=${SOURCES[1]}
             zdel=${TARGETS[0]}
             xdel=${TARGETS[1]}
             ydel=${TARGETS[2]}      
             ''')    
     
#         print nuu[j]
         SepInterp3d.PSVK3('dzK'+ref,'dxK'+ref,'dyK'+ref,
                           nuu[j],alpha[j],ref,par)
         SepInterp3d.GetPSVInK3('P'+ref,'SV'+ref,
			       'we3d-0','we3d-1','we3d-2',
                                'dzK'+ref,'dxK'+ref,'dyK'+ref,
                                ref,par)
         SepInterp3d.GetSHInK3('SH'+ref,
                               'we3d-0','we3d-1','we3d-2',
                               'dz0','dx0','dy0',
                               nuu[j],alpha[j],ref,par)
         Flow('PP'+ref,['wt'+wt,'P'+ref],
              'window n4=1 f4=%d|add ${SOURCES[1]} mode=p'%j )
         Flow('SS'+ref,['wt'+wt,'SV'+ref],
              'window n4=1 f4=%d|add ${SOURCES[1]} mode=p'%j )
         Flow('HH'+ref,['wt'+wt,'SH'+ref],
              'window n4=1 f4=%d|add ${SOURCES[1]} mode=p'%j )

    Flow('PPT-'+wt, map(lambda x: 'PP'+wt+'-%02d' % x, range(n)),
         'add ${SOURCES[1:%d]}'%n)    
    Flow('SVT-'+wt, map(lambda x: 'SS'+wt+'-%02d' % x, range(n)), 
         'add ${SOURCES[1:%d]}'%n)    
    Flow('SHT-'+wt, map(lambda x: 'HH'+wt+'-%02d' % x, range(n)), 
         'add ${SOURCES[1:%d]}'%n)  

