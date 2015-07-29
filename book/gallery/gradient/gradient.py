from rsf.proj import *
import rsf.gallery

method = rsf.gallery.method()

par = dict(xmin=2.5,xmax=7.5,zmin=0,zmax=5,
           v0=1.5,gradx=0.36,gradz=0.36,
           dim1 = 'd1=0.001 o1=0 n1=10001',
           dim2 = 'd2=0.01 o2=0 n2=501')

def igrey(custom='',title=''):
    return '''
    window min2=%g max2=%g min1=%g max1=%g |
    grey labelrot=n pclip=100 title="%s" wantaxis=y
    scalebar=n barreverse=y
    grid=y gridcol=7 screenratio=1
    label1=z unit1=km label2=x unit2=km %s
    ''' % (par['xmin'],par['xmax'],par['zmin'],par['zmax'],
           title,custom)

def zo_image(image):
    Result(image,igrey('gridcol=5','Zero-Offset %s' % method))

layers = (
    ((0,2),(3.5,2),(4.5,2.5),(5.,2.25),(5.5,2),(6.5,2.5),(10,2.5)),
    ((0,2.5),(10,3.5)),
    ((0,3.2),(3.5,3.2),(5,3.7),(6.5,4.2),(10,4.2)),
    ((0,4.5),(10,4.5))
    )

nlays = len(layers)
for i in range(nlays):
    inp = 'inp%d' % (i+1)
    Flow(inp+'.asc',None,
         'echo %s in=$TARGET data_format=ascii_float n1=2 n2=%d' % \
         (string.join(map(lambda x: string.join(map(str,x)),layers[i]),
                      ' '),
          len(layers[i])))

Flow('lay1','inp1.asc','dd form=native | spline %(dim1)s fp=0,0' % par)
Flow('lay2',None,'math %(dim1)s output="2.5+x1*0.1" ' % par)
Flow('lay3','inp3.asc','dd form=native | spline %(dim1)s fp=0,0' % par)
Flow('lay4',None,'math %(dim1)s output=4.5' % par)

Flow('lays','lay1 lay2 lay3 lay4','cat axis=2 ${SOURCES[1:4]}')
graph = '''
graph min1=%(xmin)g max1=%(xmax)g min2=%(zmin)g max2=%(zmax)g
yreverse=y wantaxis=n wanttitle=n wantscalebar=n screenratio=1
''' % par

Plot('lays0','lays',graph + ' plotfat=10 plotcol=0')
Plot('lays1','lays',graph + ' plotfat=2 plotcol=7')

# velocity
def get_velocity(vel):
    Flow(vel,None,'math %(dim1)s %(dim2)s output="%(v0)g+%(gradx)g*x1+%(gradz)g*x2" | transp' % par)
    Plot(vel,igrey('color=j allpos=y bias=1.5 title="" barlabel="v(km/s)"',par))
    Plot(vel+'-model',[vel,'lays0','lays1'],'Overlay')

Flow('dips','lays','deriv scale=y')

def zero_offset(data):
    Flow(data,'lays dips',
         '''
         kirmod twod=y freq=15 dip=${SOURCES[1]}
         dt=0.002 nt=2001
         s0=1.5 ds=0.02 ns=351
         h0=0   dh=0.02 nh=1
         type=v vel=%(v0)g gradx=%(gradx)g gradz=%(gradz)g |
         window | costaper nw2=10 | put label2=Distance unit2=km
         ''' % par)

def shots(data):
    Flow(data,'lays dips',
         '''
         kirmod twod=y freq=15 dip=${SOURCES[1]}
         dt=0.002 nt=4501
         s0=1.5 ds=0.02 ns=351
         h0=0   dh=0.04 nh=101
         type=v vel=%(v0)g gradx=%(gradx)g gradz=%(gradz)g |
         put d2=0.02 label2=Half-Offset unit2=km label3=Shot unit3=km
         ''' % par,split=[1,10001],reduce='add')

def cmps(data):
    Flow(data,'lays dips',
         '''
         kirmod twod=y freq=15 dip=${SOURCES[1]}
         dt=0.002 nt=4501 cmp=y
         s0=1.5 ds=0.02 ns=351
         h0=0   dh=0.04 nh=101
         type=v vel=%(v0)g gradx=%(gradx)g gradz=%(gradz)g |
         put d2=0.02 label2=Half-Offset unit2=km label3=Midpoint unit3=km
         ''' % par,split=[1,10001],reduce='add')

def get_impulse(impulse,data):
    Flow(impulse,data,'spike k1=751 k2=176 | smooth rect1=2 rect2=2 repeat=2')

def impulse_response(image,vel):
    Plot(image,igrey('gridcol=7','%s Impulse Response' % method))
    Plot(image+'-theory',vel,
         '''
         scale dscale=0.5 | eikonal yshot=5 | 
         contour nc=1 c0=1.5 screenratio=1 dash=1
         wantaxis=n wanttitle=n plotcol=3 plotfat=3
         min2=%g max2=%g min1=%g max1=%g
         ''' % (par['xmin'],par['xmax'],par['zmin'],par['zmax']))
    Result(image,[image,image+'-theory'],'Overlay')
