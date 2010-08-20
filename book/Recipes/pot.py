from rsf.proj import *
import fdmod,pplot,fdd,spk




def cliptogetherK(plot,file1,file2,title1,title2,axis,custom,par):

    Flow(plot+'-all',[file1,file2],
         '''
         cat axis=3 space=n ${SOURCES[1]} |
         byte gainpanel=a pclip=100|
         put label1=%(lkz)s label2=%(lkx)s
             unit1=%(ukz)s unit2=%(lkx)s
             o1=%(okz)f d1=%(dkz)f o2=%(okx)f d2=%(dkx)f 
         ''' %par)
    

    
    if(axis==1):
        if(not par.has_key('ys')): par['ys']=0.75
        if(not par.has_key('xs')): par['xs']=0.75
        if(not par.has_key('xc')): par['xc']=-8.25
        Plot(file1,plot+'-all','window n3=1 f3=0 |' + fdmod.cgrey(custom+' title=%s'%title1,par))
        Plot(file2,plot+'-all','window n3=1 f3=1 |' + fdmod.cgrey(custom+' title=%s label1=  unit1= '%title2,par))
        pplot.p1x2(plot,file1,file2,par['ys'],par['xs'],par['xc'])
    else:
        if(not par.has_key('ys')): par['ys']=.75
        if(not par.has_key('xs')): par['xs']=.75
        if(not par.has_key('xc')): par['xc']=-10
        Plot(file1,plot+'-all','window n3=1 f3=1 |' + fdmod.cgrey(custom+' title=%s wantaxis2=n label2= unit2='%title1,par))
        Plot(file2,plot+'-all','window n3=1 f3=0 |' + fdmod.cgrey(custom+' title=%s '%title2,par))
        p2x1(plot,file1,file2,par['ys'],par['xs'],par['xc'])

# ------------------------------------------------------------def cliptogetherK(plot,file1,file2,title1,title2,axis,custom,par):
def cliptogetherX(plot,file1,file2,title1,title2,axis,custom,par):
    Flow(plot+'-all',[file1,file2],
         '''
         cat axis=3 space=n ${SOURCES[1]} |
         byte gainpanel=a pclip=100 %s|
         put label1=sample label2=sample
             unit1= unit2=
                        
         ''' %(custom))
    

    
    if(axis==1):
        if(not par.has_key('ys')): par['ys']=0.75
        if(not par.has_key('xs')): par['xs']=0.75
        if(not par.has_key('xc')): par['xc']=-8.25
        Plot(file1,plot+'-all','window n3=1 f3=0 |' + fdmod.cgrey(custom+' title=%s'%title1,par))
        Plot(file2,plot+'-all','window n3=1 f3=1 |' + fdmod.cgrey(custom+' title=%s label1= unit1= '%title2,par))
        pplot.p1x2(plot,file1,file2,par['ys'],par['xs'],par['xc'])
    else:
        if(not par.has_key('ys')): par['ys']=.75
        if(not par.has_key('xs')): par['xs']=.75
        if(not par.has_key('xc')): par['xc']=-10
        Plot(file1,plot+'-all','window n3=1 f3=1 |' + fdmod.cgrey(custom+' title=%s wantaxis2=n label2= unit2='%title1,par))
        Plot(file2,plot+'-all','window n3=1 f3=0 |' + fdmod.cgrey(custom+' title=%s '%title2,par))
        p2x1(plot,file1,file2,par['ys'],par['xs'],par['xc'])

# ------------------------------------------------------------
def cliptogether(plot,file1,file2,title1,title2,axis,custom,par):

    Flow(plot+'-all',[file1,file2],
         '''
         cat axis=3 space=n ${SOURCES[1]} |
         byte gainpanel=a %s
         ''' %custom)
    

    
    if(axis==1):
        if(not par.has_key('ys')): par['ys']=0.75
        if(not par.has_key('xs')): par['xs']=0.75
        if(not par.has_key('xc')): par['xc']=-8.25
        Plot(file1,plot+'-all','window n3=1 f3=0 |' + fdmod.cgrey(custom+' title=%s'%title1,par))
        Plot(file2,plot+'-all','window n3=1 f3=1 |' + fdmod.cgrey(custom+' title=%s wantaxis1=n label1= unit1= '%title2,par))
        pplot.p1x2(plot,file1,file2,par['ys'],par['xs'],par['xc'])
    else:
        if(not par.has_key('ys')): par['ys']=.75
        if(not par.has_key('xs')): par['xs']=.75
        if(not par.has_key('yc')): par['yc']=-5
        Plot(file1,plot+'-all','window n3=1 f3=1 |' + fdmod.cgrey(custom+' title=%s wantaxis2=n label2= unit2='%title1,par))
        Plot(file2,plot+'-all','window n3=1 f3=0 |' + fdmod.cgrey(custom+' title=%s '%title2,par))
        p2x1(plot,file1,file2,par['ys'],par['xs'],par['yc']) 

def p2x1(plot,p0,p1,ys,xs,yc):
    j0 = plot + '_' + p0
    j1 = plot + '_' + p1

    Plot(j0,p0,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,0))
    Plot(j1,p1,'Overlay',vppen='yscale=%f xscale=%f ycenter=%f xcenter=0'% (ys,xs,yc))

    Plot  (plot,[j0,j1],'Overlay')
    Result(plot,[j0,j1],'Overlay')     
# ------------------------------------------------------------


def cgrey3(custom,par):
    return '''
    window min1=%g max1=%g min2=%g max2=%g min3=%g max3=%g |
    grey3 title="" framelabel=n
    label1=%s unit1=%s
    label2=%s unit2=%s
    label3=%s unit3=%s
    frame1=%d frame2=%d frame3=%d
    flat=y screenratio=%g screenht=%g point1=%g point2=%g
    %s
    ''' % (par['zmin'],par['zmax'],
           par['xmin'],par['xmax'],
           par['ymin'],par['ymax'],
           par['lz'],par['uz'],
           par['lx'],par['ux'],
           par['ly'],par['uy'],
           par['nz']/2,par['nx']/2,par['ny']/2,
           par['ratio3d'],par['height3d'],par['pointz'],par['pointx'],
           par['labelattr']+' '+custom)


def cliptogether3(plot,file1,file2,file3,title1,title2,title3,axis,center,custom,par):

    Flow(plot+'-all',[file1,file2,file3],
         '''
         cat axis=4 space=n ${SOURCES[1:3]} |
         byte gainpanel=a pclip=99.0
                 '''+custom )
    
   
    if(axis==1):
        if(not par.has_key('ys')): par['ys']=0.48
        if(not par.has_key('xs')): par['xs']=0.48
        if(not par.has_key('xc')): par['xc']=-8.5
        Plot(file1,plot+'-all','window n4=1 f4=0 | '+ cgrey3('  flat=y title=%s '%title1 +custom+center,par))
        Plot(file2,plot+'-all','window n4=1 f4=1 |' + cgrey3('  flat=y title=%s wantaxis1=n'%title2 +custom+center,par))
        Plot(file3,plot+'-all','window n4=1 f4=2 |' + cgrey3('  flat=y title=%s  wantaxis1=n'%title3 +custom+center,par))
        pplot.p1x3(plot,file1,file2,file3,par['ys'],par['xs'],par['xc'])
        Result(file1,plot+'-all','window n4=1 f4=0 | '+ cgrey3('  flat=y title=%s '%title1 +custom+center,par))
        Result(file2,plot+'-all','window n4=1 f4=1 |' + cgrey3('  flat=y title=%s '%title2 +custom+center,par))
        Result(file3,plot+'-all','window n4=1 f4=2 |' + cgrey3('  flat=y title=%s '%title3 +custom+center,par))
        
    else:
        if(not par.has_key('ys')): par['ys']=.75
        if(not par.has_key('xs')): par['xs']=.75
        if(not par.has_key('xc')): par['xc']=-10
        Plot(file1,plot+'-all','window n4=1 f4=0 | '+ cgrey3('  flat=y title=%s '%title1 +custom+center,par))
        Plot(file2,plot+'-all','window n4=1 f4=1 |' + cgrey3('  flat=y title=%s '%title2 +custom+center,par))
        Plot(file3,plot+'-all','window n4=1 f4=2 |' + cgrey3('  flat=y title=%s '%title3 +custom+center,par))
        pplot.p3x1(plot,file1,file2,file3,par['ys'],par['xs'],par['xc'])
        
  
# ------------------------------------------------------------
# extract and center wavefields
# to plot wavefield modeled with staggered grids
def displacements(uu,wfld,uz,ux,nframe,custom,par):

    nb = max(par['nb'],4)

    Flow(uz+'-t',wfld,'window n1=%d f1=%d n2=%d f2=%d n3=1 f3=0 n4=1 f4=%d'
         % (par['nz'],nb,par['nx'],nb,nframe))
    Flow(ux+'-t',wfld,'window n1=%d f1=%d n2=%d f2=%d n3=1 f3=1 n4=1 f4=%d'
         % (par['nz'],nb,par['nx'],nb,nframe))
    
    Flow(uz,uz+'-t','shift del1=%g'         % (-par['dz']/2))
    Flow(ux,ux+'-t','shift         del2=%g' % (-par['dx']/2))



# to plot wavefield modeled with centered grids
def displacementsC(uu,wfld,uz,ux,nframe,custom,par):

    nb = max(par['nb'],4)

    Flow(uz,wfld,'window n1=%d f1=%d n2=%d f2=%d n3=1 f3=0 n4=1 f4=%d'
         % (par['nz'],nb,par['nx'],nb,nframe))
    Flow(ux,wfld,'window n1=%d f1=%d n2=%d f2=%d n3=1 f3=1 n4=1 f4=%d'
         % (par['nz'],nb,par['nx'],nb,nframe))    


# ------------------------------------------------------------
# compute potentials
def potentials(rr,uz,ux,cdz,cdx,stat,custom,prefix,par):

    fdd.ddapply(rr+'zdz',uz,cdz,stat,par)
    fdd.ddapply(rr+'xdx',ux,cdx,stat,par)
    fdd.ddapply(rr+'zdx',uz,cdx,stat,par)
    fdd.ddapply(rr+'xdz',ux,cdz,stat,par)
    
    Flow(rr+'p',[rr+'zdz',rr+'xdx'],'add scale=1,+1 ${SOURCES[1]}')
    Flow(rr+'s',[rr+'xdz',rr+'zdx'],'add scale=1,-1 ${SOURCES[1]}')

#    cliptogether(rr,rr+'p',rr+'s',prefix+"P",prefix+"S",custom,par)



# compute potentials
def potentials3d(rr,uz,ux,uy,cdz,cdx,cdy,stat,custom,prefix,par):

    fdd.ddapply3d(rr+'zdz',uz,cdz,stat,par)
    fdd.ddapply3d(rr+'xdx',ux,cdx,stat,par)
    fdd.ddapply3d(rr+'ydy',uy,cdy,stat,par)
#    fdd.ddapply(rr+'zdx',uz,cdx,stat,par)
#    fdd.ddapply(rr+'xdz',ux,cdz,stat,par)
    
    Flow(rr+'p',[rr+'zdz',rr+'xdx',rr+'ydy'],'add scale=1,1,1 ${SOURCES[1:3]}')
#    Flow(rr+'s',[rr+'xdz',rr+'zdx'],'add scale=1,-1 ${SOURCES[1]}')




def plotop3d(custom,par):
    return '''
    byte gainpanel=a pclip=100 %s |
    grey3 title="" framelabel=n
    label1=%s unit1=%s
    label2=%s unit2=%s
    label3=%s unit3=%s
    flat=y screenratio=%g screenht=%g
    %s
    ''' % (custom,
           par['lz'],par['uz'],
           par['lx'],par['ux'],
           par['ly'],par['uy'],
           par['ratio3d'],par['height3d'],
           par['labelattr']+' '+custom)
