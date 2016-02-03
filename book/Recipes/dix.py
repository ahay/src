from rsf.proj import *

def dix(data,        # data name
        agc,         # migrated image
        npk,         # migration velocity
        smb,         # semblance slice
        v0,          # minimum velocity
        vm,          # median velocity
        nx,          # lateral dimension
        dx=0.01,     # midpoint sampling
        x0=0,        # lateral origin
        units='km',  # lateral units
        rect1=10,    # vertical smoothing
        rect2=10,    # lateral  smoothing
        jw=2,        # windowing for wiggle traces
        ):   
    '''Dix inversion'''

    Plot(npk,
         '''
         grey pclip=100 color=j bias=%g 
         scalebar=y title="Picked Migration Velocity"
         label1=Time unit1=s label2="Lateral Position" unit2=%s
         barlabel=Velocity barunit="%s/s" barreverse=y
         ''' % (vm,units,units))
    for line in (0,1):
        Plot(npk+str(line),npk,
             '''
             contour nc=40 wanttitle=n plotcol=%d plotfat=%d
             scalebar=y wantaxis=n barlabel=Velocity barunit="%s/s"
             ''' % ((0,10,units),(7,1,units))[line])
    Result(npk,[npk,npk+'0',npk+'1'],'Overlay')
    
    Plot(agc,
         '''
         grey title="Time-Migrated Image"
         label1=Time unit1=s label2="Lateral Position" unit2=%s
         ''' % units)
    Result(agc+'2',[agc,npk],'SideBySideAniso',vppen='txscale=1.2')
    
    slp = data+'-slp'
    Flow(slp,agc,'dip liter=40 rect1=%d rect2=%d' % (rect1,rect2))
    Plot(slp,
         '''
         grey title="Estimated Dips" scalebar=y color=j
         label1=Time unit1=s label2="Lateral Position" unit2=%s pclip=100
         barlabel="Dip (samples)"
         ''' % units)
    Result(slp,[agc,slp],'SideBySideAniso')

    vg = data+'-vg'
    Flow(vg,[agc,slp],
         'pwdsigk dips=${SOURCES[1]} verb=y nliter=5 niter=100')
    Plot(vg,
         '''
         grey title="Structure Enhancing"
         label1=Time unit1=s label2="Lateral Position" unit2=%s
         ''' % units)
    Plot(vg+'2',[agc,vg],
         '''
         add scale=1,-1 ${SOURCES[1]} | cat ${SOURCES[0]} |
         byte gainpanel=2 | window n3=1 |
         grey title="Removed Noise" 
         label1=Time unit1=s label2="Lateral Position" unit2=%s
         ''' % units)    
    Result(vg,[vg,vg+'2'],'SideBySideAniso')

    
    dix = data+'-dix'
    Flow([dix,dix+'0'],[npk,smb],
         '''
         dix rect1=%d rect2=%d
         weight=${SOURCES[1]} vrmsout=${TARGETS[1]}
         niter=30
         ''' % (rect1/2,rect2/2))
    Plot(dix,
         '''
         grey pclip=100 color=j bias=%g 
         scalebar=y title="Interval Velocity" barlabel=Velocity barunit="%s/s"
         label1=Time unit1=s label2="Lateral Position" unit2=%s barreverse=y
         ''' % (vm,units,units))
    Plot(dix+'0',
         '''grey pclip=100 color=j bias=%g 
         scalebar=y title="Predicted Migration Velocity" barlabel=Velocity barunit="%s/s"
         label1=Time unit1=s label2="Lateral Position" unit2=%s barreverse=y
         ''' % (vm,units,units))
    Result(dix,[dix,dix+'0'],'SideBySideAniso')
    Result(dix+'0',[npk,dix+'0'],'SideBySideAniso')

#    pdx = data+'-pdx'
#    Flow([pdx,pdx+'0'],[npk,smb,slp],
#         '''
#         pwdix slope=${SOURCES[2]}
#         weight=${SOURCES[1]} vrmsout=${TARGETS[1]}
#         niter=50 verb=y ncycle=10 rect1=%d
#         ''' % (4*rect1))
#    Plot(pdx,
#         '''
#         grey pclip=100 color=j bias=%g allpos=y
#         scalebar=y barlabel="Estimated Interval Velocity" barunit="%s/s"
#         label1=Time unit1=s label2="Lateral Position" unit2=%s
#         wanttitle=n 
#         ''' % (v0,units,units))
#    Plot(pdx+'0',
#         '''
#         grey pclip=100 color=j bias=%g
#         scalebar=y title="Predicted Migration Velocity"
#         label1=Time unit1=s label2="Lateral Position" unit2=%s 
#         barlabel=Velocity barunit="%s/s" barreverse=y
#         ''' % (vm,units,units))
#    Result(pdx,[npk,pdx+'0'],'SideBySideAniso')

    shp = data+'-shp'
    ext = data+'-ext'

    for inp in (npk,agc+'2',slp):
        Flow(inp+'_b',inp,'window n2=1 | spray axis=2 n=10 o=%g d=%g' % (x0-10*dx,dx))
        Flow(inp+'_e',inp,'window n2=1 f2=%d | spray axis=2 n=10' % (nx-1))
        Flow(inp+'_',[inp+'_b',inp,inp+'_e'],'cat ${SOURCES[1:3]} axis=2')
    
    Flow([shp,shp+'0'],[npk+'_',agc+'2_',slp+'_'],
         '''
         dixshape dip=${SOURCES[2]}
         weight=${SOURCES[1]} vrmsout=${TARGETS[1]}
         niter=50 verb=y rect1=%d rect2=5 lam=0.1
         ''' % rect1)
    Plot(shp,
         '''
         window f2=10 n2=%d |
         grey pclip=100 color=j bias=%g allpos=y
         scalebar=y barlabel="Estimated Interval Velocity" barunit="%s/s"
         label1=Time unit1=s label2="Lateral Position" unit2=%s barreverse=y
         wanttitle=n 
         ''' % (nx,v0,units,units))
    Plot(shp+'1',shp,
         '''
         window f2=10 n2=%d |
         grey pclip=100 color=j bias=%g allpos=y
         scalebar=y barlabel="Velocity" barunit="%s/s"
         label1=Time unit1=s label2="Lateral Position" unit2=%s barreverse=y
         title="Interval Velocity" 
         ''' % (nx,v0,units,units))
    Plot(shp+'0',
         '''
         window f2=10 n2=%d |
         grey pclip=100 color=j bias=%g
         scalebar=y title="Predicted MIgration Velocity"
         label1=Time unit1=s label2="Lateral Position" unit2=%s 
         barlabel=Velocity barunit="%s/s" barreverse=y
         ''' % (nx,vm,units,units))
    Result(shp,[shp+'1',shp+'0'],'SideBySideAniso')
    
    Plot(agc+'w',agc,
         '''
         window j2=%d |
         wiggle transp=y yreverse=y scalebar=y wantaxis=n wanttitle=n
         plotcol=7 poly=y
         ''' % jw)
    
    Plot(vg+'w',vg,
         '''
         window j2=%d |
         wiggle transp=y yreverse=y scalebar=y wantaxis=n wanttitle=n
         plotcol=7 poly=y
         ''' % jw)
    
    Result(dix+'w',[dix,agc+'w'],'Overlay')
#    Result(pdx+'w',[pdx,vg+'w'],'Overlay')
    Result(shp+'w',[shp,agc+'w'],'Overlay')
