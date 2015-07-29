from rsf.proj import *

# ================================= #
# model parameters                  #
# ================================= #

ft = .3048

nz = 1201; z0 = 0.; dz = 25*ft
nx = 3201; x0 = 0.; dx = 25*ft
nT = 751;  T0 = 0.; dT = .0067
nt = 37501;t0 = 0.; dt = .00032

ns = 500;  s0 = 925 * ft; ds = 150*ft; zs = 25*ft; Ts = 5.08e-3
nh = 348;  h0 = 0.;       dh = 75*ft;  zr = 25*ft; Tr = 5.08e-3

nC = 3;    C0 = 800;      dC = 800  # cig position
nH = 64;                  dH = 1    # subsurface offset

z = map(lambda t: z0 + t * dz, range(nz))
x = map(lambda t: x0 + t * dx, range(nx))
T = map(lambda t :T0 + t * dT, range(nT))
t = map(lambda t: t0 + t * dt, range(nt))
s = map(lambda t: s0 + t * ds, range(ns))
h = map(lambda t: h0 + t * dh, range(nh))

freq = 20
j3 = 15; n3 = 1 + int((nt-1)/j3); d3 = dt * j3
j4 = 25; n4 = 1 + int((nt-1)/j4); d4 = dt * j4

b = [40   ,40   ,40   ,40   ]
c = [4.e-5,2.e-5,2.e-5,2.e-5]
eps_sigm = 1

shot0 = 0
nshot = 1
jshot = 1
ss = []

job_name = ''
job_path = ''

bzl = b[0]; bzh = b[1]; bxl = b[2]; bxh = b[3]; bTl = 0; bTh = 0
czl = c[0]; czh = c[1]; cxl = c[2]; cxh = c[3]; cTl = 0; cTh = 0

Nx = nx + bxl + bxh
Nz = nz + bzl + bzh
NT = nT + bTl + bTh

# half opening angle for cig
# need abs(p[:]) >= tan(a[-1])
na = 151; a0 = 0;    da = .4 # deg
np = 401; p0 = -1.8; dp = .009

# ================================= #
# preprocessing                     #
# ================================= #

def mapping(out,inp,tau,inv=False):
    '''z to tau : inv=False
    tau to z : inv=True'''
    if inv:
        Flow(out,[inp,tau],'''
        pseudodepth tau=${SOURCES[1]}
        inv=y linear=y 
        ''')
    else:
        Flow(out,[inp,tau],'''
        pseudodepth tau=${SOURCES[1]}
        inv=n n=%d o=0 d=%g linear=y
        ''' % (nT,dT))

# velocity
def get_velo(velo):
    Fetch('sigsbee2a_migration_velocity.sgy','sigsbee')
    Flow([velo,velo+'_t'],'sigsbee2a_migration_velocity.sgy','''
    segyread tfile=${TARGETS[1]} hfile=/dev/null bfile=/dev/null |
    put o1=0 d1=25 o2=10025 d2=37.5 | transp |
    remap1 n1=3201 o1=10000 d1=25 | transp | scale rscale=.3048 |
    put o1=%g d1=%g label1=z unit1=m
        o2=%g d2=%g label2=x unit2=m |
    clip2 lower=1450
    ''' % (z0,dz,x0,dx))

def grey_velo(velo):
    Result(velo,'grey scalebar=y barreverse=y bias=1450 color=j title=velocity')

def grey_velo4(velo,Velo,VElo):
    vgrey = 'grey scalebar=y bias=1500 color=j allpos=y'
    Plot(velo,vgrey + ' title=car')
    Plot(Velo,vgrey + ' title=tau')
    Plot(VElo,vgrey + ' title=back')
    Plot('verr',[velo,VElo],'add scale=-1,1 ${SOURCES[1]} | grey color=j scalebar=y title=error')
    Result('verr',[velo,Velo,VElo,'verr'],'TwoRows')

def get_vmap(vmap,velo,smooth='smooth rect1=16 rect2=16 repeat=8'):
    '''mapping velociy'''
    if not len(smooth):
        smooth = 'math output=input'
    Flow(vmap,velo,smooth)

def get_tau(tau,vmap):
    '''vertical time'''
    Flow(tau,vmap,'math output="1./input" | integral1 rule=s')

def cont_tau(tau,velo):
    '''velo can be sharp or smoothed'''
    Plot(tau,'contour c0=0 dc=.4 allpos=y plotfat=10 title="tau (dc=0.4sec)"')
    Plot('velo1',velo,'grey allpos=y bias=1500 wantaxis=n wanttitle=n')
    Result(tau,['velo1',tau],'Overlay')

def get_sigm(sigm,tau):
    Flow(sigm,tau,'transp | deriv scale=y | transp')

def grey_sigm2(sigm,Sigm):
    Plot(sigm,'grey color=j scalebar=y')
    Plot(Sigm,'grey color=j scalebar=y')
    Result(sigm,[sigm,Sigm],'SideBySideIso')

def get_sour(sour):
    Flow(sour,None,'''
    spike nsp=1 mag=1 n1=%d o1=0 d1=%g k1=%d | 
    ricker1 frequency=%g | scale axis=1 | put label1=t unit1=sec
    ''' % (nt,dt,1+1.0/(freq*dt),freq))
def graph_sour(sour,n=800):
    Result(sour,'''window n1=%d |
    graph title="%g Hz, total %g sec" plotfat=10
    ''' % (n,freq,t[-1]))

def set_shot(shot0_,nshot_,jshot_):
    global shot0,nshot,jshot,ss
    shot0 = shot0_
    nshot = nshot_
    jshot = jshot_
    ss = s[shot0:shot0+nshot*jshot:jshot]

# full dataset
def get_dataset(dataset):
    Fetch('sigsbee2a_nfs.sgy','sigsbee')
    Flow('sigsbee2a_nfs hdr','sigsbee2a_nfs.sgy','''
    segyread tfile=hdr hfile=/dev/null bfile=/dev/null > ${TARGETS[0]} &&
    dd < hdr type=float > ${TARGETS[1]} && rm hdr
    ''',stdout=-1)
    Flow('hdr_ss','hdr','headermath output="10925/150+fldr" | window')
    Flow('hdr_rr','hdr','headermath output="offset/75"      | window')
    Flow('hdr_ssrr','hdr_rr hdr_ss','cat ${SOURCES[1]} axis=2 | transp | dd type=int')
    Flow(dataset,'sigsbee2a_nfs hdr_ssrr','''
    intbin head=${SOURCES[1]} xkey=0 ykey=1 |
    put             label1=t unit1=sec
        o2=%g d2=%g label2=offset unit2=m
        o3=%g d3=%g label3=shot unit3=m
    ''' % (h0,dh,s0,ds))

# shot gathers in current job
def get_data(data,dataset):
    Flow(data,dataset,'''
    pad end1=1 |
    window n3=%d f3=%d j3=%d n1=%d |
    remap1 n1=%d o1=0 d1=%g
    ''' % (nshot,shot0,jshot,n4,nt,dt))

def grey_data(data):
    Result(data,'''
    window j1=%d | byte | transp plane=23 memsize=16384 |
    grey3 flat=n frame1=250 frame3=0 frame2=64 point1=.6 point2=.6
    ''' % j4)

def plot_survey(velo):
    limit = ' min1=%g max1=%g min2=%g max2=%g' % (x[0],x[-1],z[0],z[-1])
    
    Plot('vb',velo,'grey color=I wanttitle=n')
    for i in range(0,len(ss),50):
        cshot = ss[i] # current shot
        hmin = h[0]   # current min offset
        hmax = h[-1]  # current max offset
        if cshot + hmin > x[-1]:
            continue
        if cshot + hmax > x[-1]:
            nchan = int((x[-1] - cshot - h[0]) / dh) + 1
            hmax = (nchan-1) * dh + h[0]
        else:
            nchan = nh
        # plot source
        Plot('ss%d ss%d.asc' % (i,i),None,'''
        echo "%g %g" > ${TARGETS[1]} &&
        echo "n1=2 data_format=ascii_float in=${TARGETS[1]}" |
        dd type=complex |
        graph title="shot %d x=%g z=%g nchan=%d" wheretitle=b
        symbol=* symbolsz=15 plotcol=5 plotfat=10
        wantaxis=n yreverse=y
        ''' % (cshot,zs,i,cshot,zs,nchan) + limit)
        # plot receiver line
        Flow('x%d x%d.asc' % (i,i),None,'''
        for j in {%d..%d}; do echo %g+%g*$$j | bc; done > ${TARGETS[1]} &&
        echo "n1=%d in=${TARGETS[1]} data_format=ascii_float" |
        dd form=native
        ''' % (0,nchan-1,cshot+hmin,dh,nchan))
        Flow('z%d z%d.asc' % (i,i),None,'''
        for j in {%d..%d}; do echo %g; done > ${TARGETS[1]} &&
        echo "n1=%d in=${TARGETS[1]} data_format=ascii_float" |
        dd form=native
        ''' % (0,nchan-1,zs,nchan))
        Plot('rr%d' % i,'x%d z%d' % (i,i),'''
        cmplx ${SOURCES[1]} |
        graph plotfat=16 plotcol=6 wantaxis=n wanttitle=n yreverse=y
        ''' + limit)
        Plot('ssrr%d' % i,'vb rr%d ss%d' % (i,i),'Overlay')
        Plot('shot%d' % i,'data','window n3=1 f3=%d j1=%d | grey ' % (i,j4))
        Result('ssrr%d' % i,['ssrr%d' %i,'shot%d' % i],'SideBySideIso')

# ================================= #
# job scripts                       #
# ================================= #

def set_job(job='test'):
    global job_name,job_path
    job_name = job;
    job_path = os.path.join(os.environ['DATAPATH'],job)
    if not os.path.isdir(job_path):
        os.mkdir(job_path)

def swapbyte(files):
    for source in files:
        target = source + '_job'
        fname = os.path.join(job_path,source.lower() + '.bin')
        Flow(target,source,'swapbyte verb=y out=%s' %  fname)

def make_par(forward=True,backward=True,pseudo=False):
    mode = 0
    if backward:
        mode += 2
    if forward:
        mode += 1
    mode = mode % 3
    # savedata=False hard-coded

    f = open(os.path.join(job_path,'par.txt'),'w')
    if pseudo:
        f.write('tau=1\nnz=%d\nz0=%g\ndz=%g\nzs=%g\nzr=%g\neps=%g\n' % (nT,T0,dT,Ts,Tr,eps_sigm))
    else:
        f.write('tau=0\nnz=%d\nz0=%g\ndz=%g\nzs=%g\nzr=%g\n' % (nz,z0,dz,zs,zr))
    f.write('''mode=%d
nx=%d\nx0=%g\ndx=%g
nt=%d\ndt=%g\nj3=%d
bzl=%d\nbzh=%d\nbxl=%d\nbxh=%d
czl=%g\nczh=%g\ncxl=%g\ncxh=%g
ns=%d\ns0=%g\nds=%g
nh=%d\nh0=%g\ndh=%g
nH=%d\ndH=%d
nC=%d\ndC=%d\nC0=%d
''' % (mode,nx,x0,dx,nt,dt,j3,
           b[0],b[1],b[2],b[3],
           c[0],c[1],c[2],c[3],
           nshot,s[shot0],ds*jshot,
           nh,   h0,      dh,
           nH,dH,nC,dC,C0))
    f.close()

def make_job(nproc=64,nnode=64,pseudo=False):
    bg_mode = ['SMP','DUAL','VN'][int(nproc/(2*nnode))] # should be log2()

    f = open(os.path.join(job_path,'job.sh'),'w+')
    f.write('''#!/bin/sh
# @ account_no = k09
# @ job_name = %s
# @ job_type = bluegene
# @ output = $(job_name).out
# @ error = $(job_name).err
# @ environment = COPY_ALL
# @ wall_clock_limit = 23:59:00
# @ notification = never
# @ bg_size = %d
# @ queue

cd %s
mpirun -mode %s -np %d -cwd $PWD -exe ./rtmiso.exe
''' % (job_name,nnode,job_path,bg_mode,nproc))
    f.close()

# ================================= #
# postprocessing                    #
# ================================= #

# update bzl bxl etc. fpar is par.txt
def update_par(fpar,pseudo):
    global bzl,bzh,bxl,bxh,bTl,bTh,Nx,Nz,NT,ns,ns,d3

    f = open(fpar,'r')
    lines = f.readlines()
    f.close()

    par = {}
    for line in lines:
        key,val = line.replace('out : ','')[:-1].split('=')
        par[key] = val
    Nx,bxl,bxh = map(lambda t: int(par[t]), ['Nx','bxl','bxh'])
    if pseudo:
        NT,bTl,bTh = map(lambda t: int(par[t]), ['Nz','bzl','bzh'])
    else:
        Nz,bzl,bzh = map(lambda t: int(par[t]), ['Nz','bzl','bzh'])
    n3 = int(par['n3'])
    d3 = float(par['d3'])
    ns = int(par['ns'])

# wavefield movie
def grey_wave(i,sr='s',pseudo=False):
    '''sr={s,r}'''
    endian = ' > tmp.rsf && swapbyte < tmp.rsf verb=y'

    if pseudo:
        target = 'Uu' + sr + '%d' % i
        
        Flow(target,None,'''
        echo "n1=%d o1=%g d1=%g label1=tau unit1=sec
              n2=%d o2=%g d2=%g label2=x unit2=m
              n3=%d o3=0  d3=%g label3=t unit3=sec
              in=%s.bin data_format=native_float"
        ''' % (NT,-bTl*dT,dT,
               Nx,-bxl*dx,dx,
               n3,d3,target+'BE') + endian)
    else:
        target = 'uu' + sr + '%d' % i

        Flow(target,None,'''
        echo "n1=%d o1=%g d1=%g label1=z unit1=m
              n2=%d o2=%g d2=%g label2=x unit2=m
              n3=%d o3=0  d3=%g label3=t unit3=sec
              in=%s.bin data_format=native_float"
        ''' % (Nz,-bzl*dz,dz,
               Nx,-bxl*dx,dx,
               n3,d3,target+'BE') + endian)
    Result(target,'window j3=8 | grey gainpanel=a')

# stack partial images
def make_imag(imag,pseudo=False):
    endian = ' > tmp.rsf && swapbyte < tmp.rsf verb=y'

    imags = []
    for i in range(ns):
        if pseudo:
            target = 'Img%d' % i

            Flow(target,None,'''
            echo "n1=%d o1=%g d1=%g label1=tau unit1=sec
                  n2=%d o2=%g d2=%g label2=x unit2=m
                  in=%s.bin data_format=native_float"
            ''' % (NT,-bTl*dT,dT,
                   Nx,-bxl*dx,dx,
                   target+'BE') + endian)
        else:
            target = 'img%d' % i

            Flow(target,None,'''
            echo "n1=%d o1=%g d1=%g label1=z unit1=m
                  n2=%d o2=%g d2=%g label2=x unit2=m
                  in=%s.bin data_format=native_float"
            ''' % (Nz,-bzl*dz,dz,
                   Nx,-bxl*dx,dx,
                   target+'BE') + endian)
        imags.append(target)

    if pseudo:
        cut = ' | window n1=%d f1=%d n2=%d f2=%d' % (nT,bTl,nx,bxl)
    else:
        cut = ' | window n1=%d f1=%d n2=%d f2=%d' % (nz,bzl,nx,bxl)
    Flow(imag,imags,'cat axis=3 ${SOURCES[1:%d]} | stack axis=3' % len(imags) + cut + ' out=%s.bin' % imag)

# stack partial cigs
def make_cig(cig,pseudo=False):
    endian = ' > tmp.rsf && swapbyte < tmp.rsf verb=y'

    cigs = []
    for i in range(ns):
        if pseudo:
            target = 'Cig%d' % i

            Flow(target,None,'''
            echo "n1=%d o1=%g d1=%g label1=tau unit1=sec
                  n2=%d o2=%g d2=%g label2=x unit2=m
                  n3=%d o3=%g d3=%g label3=h unit3=m
                  in=%s.bin data_format=native_float"
            ''' % (NT,-bTl*dT,dT,
                   nC,x[C0],dx*dC,
                   nH,0,dx*dH,
                   target+'BE') + endian)
        else:
            target = 'cig%d' % i

            Flow(target,None,'''
            echo "n1=%d o1=%g d1=%g label1=z unit1=m
                  n2=%d o2=%g d2=%g label2=x unit2=m
                  n3=%d o3=%g d3=%g label3=h unit3=m
                  in=%s.bin data_format=native_float"
            ''' % (Nz,-bzl*dz,dz,
                   nC,x[C0],dx*dC,
                   nH,0,dx*dH,
                   target+'BE') + endian)
        cigs.append(target)

    if pseudo:
        cut = ' | window n1=%d f1=%d' % (nT,bTl)
    else:
        cut = ' | window n1=%d f1=%d' % (nz,bzl)
    Flow(cig,cigs,'cat axis=4 ${SOURCES[1:%d]} | stack axis=4' % len(cigs) + cut + ' out=%s.bin' % cig)

def grey_imag(imag,pseudo):
    if pseudo:
        gain = 'pow pow1=2'
        filt = 'gradtest grad=2 | bandpass flo=10.804 fhi=22.746'
        label= 'put unit2=km d2=%g' % (.001*dx)
    else:
        gain = 'pow pow1=2'
        filt = 'gradtest grad=2 | bandpass flo=.0095 fhi=.02'
        label= 'put d1=%g unit1=km d2=%g unit2=km' % (.001*dz,.001*dx)
	# note here ktau/kz = dz/dtau = 1137.3

    Result(imag,imag,gain + ' | ' + filt + ' | ' + label + ' | grey color=i pclip=98 screenwd=12 screenht=4 labelfat=2 labelsz=3 title=')

# offset to angle (cartesian only)
def make_acig(acig,cig):
    Flow(acig,cig,'''
    transp plane=23 |
    slant adj=y np=%d p0=%g dp=%g verb=y |
    tan2ang na=%d a0=%g da=%g | 
    put label2=angle unit2=deg out=%s.bin
    ''' % (np,p0,dp,na,a0,da,acig))

# plotting angle cig
def grey_acig(acig):
    for i in range(nC):
        target = acig + '%d' % i
        Result(target,acig,'''
        window n3=1 f3=%d |
        pow pow1=1 | bandpass flo=.01 |
        put d1=%g unit1=km |
        grey labelfat=4 labelsz=5 titlefat=4 titlesz=8 title="%g km" screenwd=6 screenht=11
        ''' % (i,.001*dz,.001*x[C0+i*dC]))

# plotting angle cig in one figure
def grey_acig2(acig):
	Result(acig,'put n2=%d n3=1 d1=%g unit1=km | pow pow1=1 | bandpass flo=.01 | grey labelfat=2 labelsz=3 screenwd=12 screenht=4 wantaxis2=n title= pclip=99' % (na*nC,.001*dz))

# cig(tau,x,h) to cig(z,x,h)
def cig_tau2z(CIg,Cig,tau):
    Flow('tau3',tau,'window n2=%d f2=%d d2=%d | spray axis=3 n=%d o=%d d=%d' % (nC,C0,dC,nH,0,dH))
    mapping(CIg,Cig,'tau3',inv=True)

# x axis origin has been shifted to zero
# nx = 3201; x0 = 10000*ft; dx = 25*ft
# ns = 500;  s0 = 10925*ft; ds = 150*ft; zs = 25*ft; Ts = 5.08e-3
