#
# Functions for making deep-water channel models.
#
# James W. Jennings Jr.
# Research Scientist
# 
# Bureau of Economic Geology
# John A. and Katherine G. Jackson School of Geosciences
# University of Texas at Austin
# University Station, Box X
# Austin, TX 78713-8924
# 
# 512-471-1534 (voice)
# 512-471-0140 (fax)
# mailto:jim.jennings@beg.utexas.edu
# http://www.beg.utexas.edu/staffinfo/jennings01.htm
# 
# July 2005
#
# $Id$
#

#
# Import RSF and math libraries.
#

from rsfproj import *
from math    import *
import rfield

#
# Convert props.dat format to RSF ascii
#

def format_props (target=None,source=None,env=None):

    inp   = open(str(source[0]),'r')    # open files
    file  = str(target[0])
    out   = open(file,'w')
    
    lines = inp.readlines()             # read input file
    lines.pop(0)                        # remove the first line
    for line in lines:                  # get all fields except the first one
        out.write(string.join(line.split()[1:],' ')+'\n')
        
    out.write('''
              n1=4 d1=1 o1=1
              n2=%d d2=1 o2=1
              esize=0 in=%s data_format=ascii_float
              ''' % (len(lines),file) )
              
    out.close()                         # close files
    inp.close()

#
# Get channel properties from the data server & convert to RSF
#

def get_props (private=None):

    for type in ('skew','pos'):         # get skew and pos arrays from server
        Fetch (type+'_data.rsf','jim',private)
    
    Fetch ('props.dat','jim',private)   # get props table from server
                                        # make a builder
    Props = Builder(action=Action(format_props),src_suffix='.dat',suffix='.asc')
    env   = Environment(BUILDERS={'Props':Props})
    
    env.Props('props')                  # run the builder(?)
                                        # convert ascii to native
    Flow ('props','props.asc','dd form=native')

#
# Put the channel properties into layer arrays
#

def make_prop_arrays (  memsize=100,
                        nx=None, ny=None,
                        dx=None, dy=None,
                        ox=None, oy=None):

    for type in ('skew','pos'):         # interpolate skew and pos arrays
        Flow (type,type+'_data',
              '''
              dd form=native |
              remap1 n1=%d d1=%g o1=%g | transp memsize=%d plane=12 |
              remap1 n1=%d d1=%g o1=%g | transp memsize=%d plane=12
              ''' % (nx,dx,ox,memsize,ny,dy,oy,memsize) )

                                        # get the props out of the table
    Flow ('depth_list'  ,'props','window f1=0 n1=1')
    Flow ('wtop_list'   ,'props','window f1=1 n1=1')
    Flow ('wbot_list'   ,'props','window f1=2 n1=1')
    Flow ('aggrade_list','props','window f1=3 n1=1 | add add=20')

                                        # put the props in layer arrays
    for prop in ('depth','wtop','wbot','aggrade'):
        Flow (prop,prop+'_list',
              '''
              spray axis=2 n=%d d=%g o=%g | 
              spray axis=3 n=%d d=%g o=%g |
              transp memsize=%d plane=13
              ''' % (ny,dy,oy,nx,dx,ox,memsize) )

#
# Make channel profile surfaces
#

def make_surfaces ( bd_depth=None,  md_depth=None,
                    bd_zshift=None, md_zshift=None):

    Flow ('xleft','pos skew wtop wbot',
          '''
          math p=${SOURCES[0]} s=${SOURCES[1]} wt=${SOURCES[2]} wb=${SOURCES[3]} 
          output="(2*p+s*(wt-wb)-wb)/((s+1)*(wt-wb))"
          ''',stdin=0)
        
    Flow ('xrght','pos skew wtop wbot',
          '''
          math p=${SOURCES[0]} s=${SOURCES[1]} wt=${SOURCES[2]} wb=${SOURCES[3]} 
          output="(2*p+s*(wt-wb)+wb)/((s-1)*(wt-wb))"
          ''',stdin=0)
         
    prefix = ["ch_" ,"bd_"      ,"md_"]
    depthf = [ 1    , bd_depth  , md_depth]
    zshift = [ 0    , bd_zshift , md_zshift]
    
    for i in (0,1,2):
        Flow (prefix[i]+'profiles','depth xleft xrght',
              '''
              math d=depth.rsf xl=xleft.rsf xr=xrght.rsf
              output="d*(%g+%g*(tanh(%g*(xl-0.5))+tanh(%g*(xr-0.5)))/2)" |
              clip2 upper=0
              ''' % (zshift[i],depthf[i],pi,pi) )
    
    Flow ('ch_surfaces','ch_profiles aggrade','add ${SOURCES[1]}')
    
    Flow ('bd_surfaces','bd_profiles ch_profiles aggrade',
          '''
          minmax mode=max file1=${SOURCES[0]} file2=${SOURCES[1]} |
          add ${SOURCES[2]}
          ''')
    
    Flow ('md_surfaces','md_profiles ch_profiles aggrade',
          '''
          minmax mode=max file1=${SOURCES[0]} file2=${SOURCES[1]} |
          add ${SOURCES[2]}
          ''')

#
# Erode the channel profile surfaces
#

def erode_surfaces (memsize=100):

    Flow ('ch_erode','ch_surfaces',
          '''
          transp memsize=%d plane=13 | reverse which=1 | 
          listminmax | 
          reverse which=1 | transp memsize=%d plane=13
          ''' % (memsize,memsize) )
          
    Flow ('ch_erode_temp','ch_erode','window squeeze=n f3=1 n3=14')
    
    Flow ('bd_erode','bd_surfaces ch_erode_temp',
          '''
          window squeeze=n f3=14 n3=1 > bd_surface_15.rsf &&
          cat axis=3 ${SOURCES[1]} bd_surface_15.rsf > bd_erode_temp.rsf &&
          minmax mode=min file1=bd_erode_temp.rsf file2=${SOURCES[0]} > $TARGET &&
          rm bd_surface_15.rsf bd_erode_temp.rsf
          ''', stdout=0)

    Flow ('md_erode','md_surfaces ch_erode_temp',
          '''
          window squeeze=n f3=14 n3=1 > md_surface_15.rsf &&
          cat axis=3 ${SOURCES[1]} md_surface_15.rsf > md_erode_temp.rsf &&
          minmax mode=min file1=md_erode_temp.rsf file2=${SOURCES[0]} > $TARGET &&
          rm md_surface_15.rsf md_erode_temp.rsf
          ''', stdout=0)

#
# Put channel properties into a regular grid
#

def make_prop_grids (memsize=100,nz=None,dz=None,oz=None):

    for prop in ('pos','skew','aggrade'):
        Flow (prop+'_pad',prop,'pad beg3=1 o3=0')
        Flow (prop+'_grid','ch_erode '+prop+'_pad',
              '''
              unif3 n1=%d d1=%g o1=%g layers=${SOURCES[1]} | 
              transp memsize=%d plane=13 | transp memsize=%d plane=12
              ''' % (nz,dz,oz,memsize,memsize) )
    
    for prop in ('depth','wtop','bd_erode','md_erode'):
        Flow (prop+'_pad' ,prop,
              'window squeeze=n f3=0 n3=1 | cat axis=3 $SOURCE')
        Flow (prop+'_grid','ch_erode '+prop+'_pad',
              '''
              unif3 n1=%d d1=%g o1=%g layers=${SOURCES[1]} | 
              transp memsize=%d plane=13 | transp memsize=%d plane=12
              ''' % (nz,dz,oz,memsize,memsize) )

#
# Make indicators and masks
#

def make_masks (memsize=100,nz=None,dz=None,oz=None):

# layer indicator

    Flow ('layer','ch_erode',
          '''
          unif3 n1=%d d1=%g o1=%g v00=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 |
          transp memsize=%d plane=13 | transp memsize=%d plane=12
          ''' % (nz,dz,oz,memsize,memsize) )
    
# above-layer-zero mask

    Flow ('layer_00_mask','layer','clip2 upper=1')
    
# channel mask

    Flow ('ch_mask','z_rel','mask max=-0.01 | dd type=float')

# bypass drape and margin drape masks

    if (nz > 1):
        Flow ('bd_mask','bd_erode_grid ch_mask',
              '''
              math surf=${SOURCES[0]} output="x3-surf" | 
              mask max=0 | dd type=float | add mode=m ${SOURCES[1]}
              ''')
        
        Flow ('md_mask','md_erode_grid ch_mask',
              '''
              math surf=${SOURCES[0]} output="x3-surf" | 
              mask max=0 | dd type=float | add mode=m ${SOURCES[1]}
              ''')
    else:
        Flow ('bd_mask','bd_erode_grid ch_mask',
              '''
              math surf=${SOURCES[0]} output="%g-surf" | 
              mask max=0 | dd type=float | add mode=m ${SOURCES[1]}
              ''' % (oz) )
              
        Flow ('md_mask','md_erode_grid ch_mask',
              '''
              math surf=${SOURCES[0]} output="%g-surf" | 
              mask max=0 | dd type=float | add mode=m ${SOURCES[1]}
              ''' % (oz) )

# amalgamated sand mask

    Flow ('as_mask','sand_parm bd_mask md_mask',
          '''
          mask min=1 | dd type=float |
          math bd=${SOURCES[1]} md=${SOURCES[2]} output="input*(1-bd)*(1-md)"
          ''')
    
# non-amalgamated sand mask

    Flow ('na_mask','ch_mask bd_mask md_mask as_mask',
          'add scale=1,-1,-1,-1 ${SOURCES[1]} ${SOURCES[2]} ${SOURCES[3]}')

#
# Make sand fraction grid
#

def make_sand ( memsize=100,
                nz=None, dz=None, oz=None,
                as_width=None,   as_shift=None, 
                as_height0=None, as_height1=None, as_shape=None,
                bd_sand=None,    md_sand=None,    as_sand=None, 
                na_sand0=None,   na_sand1=None):

    if (nz > 1):
        Flow ('z_rel','aggrade_grid depth_grid layer_00_mask',
              '''
              math a=${SOURCES[0]} d=${SOURCES[1]} m=${SOURCES[2]}
              output="m*(x3-a)/d" |
              clip2 upper=0
              ''')
    else:
        Flow ('z_rel','aggrade_grid depth_grid layer_00_mask',
              '''
              math a=${SOURCES[0]} d=${SOURCES[1]} m=${SOURCES[2]}
              output="m*(%g-a)/d" |
              clip2 upper=0
              ''' % (oz) )
         
    Flow ('pos_shift','pos_grid wtop_grid skew_grid',
          '''
          math p=${SOURCES[0]} w=${SOURCES[1]} s=${SOURCES[2]}
          output="p/w+%g*s/2"
          ''' % (as_shift) )
         
    Flow ('as_height','skew_grid',
          'math s=$SOURCE output="%g+%g*s"'
          % (as_height0,as_height1-as_height0) )
         
    Flow ('sand_parm','pos_shift z_rel as_height',
          '''
          math x=${SOURCES[0]} z=${SOURCES[1]} h=${SOURCES[2]}
          output="-z/(1-h+4^%g*h*(x^2/%g^2)^%g)"
          ''' % (as_shape,as_width,as_shape) )
         
    Flow ('sand_parm_clip','sand_parm','clip2 upper=1')
    
    Flow ('sand','sand_parm_clip bd_mask md_mask as_mask na_mask',
          '''
          math p=${SOURCES[0]} 
          bd=${SOURCES[1]} md=${SOURCES[2]} as=${SOURCES[3]} na=${SOURCES[4]}
          output="bd*%g+md*%g+as*%g+na*(%g+%g*p)"
          ''' % (bd_sand,md_sand,as_sand,na_sand0,na_sand1-na_sand0) )
     
#
# Make shale and porosity grids
#

def make_shale_and_phi (prefix=''):

    Flow (prefix+'shale',prefix+'sand_xypad',
          'math sand=$SOURCE output="1-sand"')
    
    Flow (prefix+'phi',prefix+'shale',
          '''
          math shale=$SOURCE
          output="0.378-0.144*shale-0.0000262*(10000-7000)"
          ''')

#
# Make density and sonic velocities
#

def make_rho_and_vel (prefix=''):

    Flow (prefix+'rho_noise',[prefix+'phi_noise',prefix+'shale'],
          '''
          math phi=${SOURCES[0]} shale=${SOURCES[1]}
          output="2.70-0.085*shale-1.75*phi"
          ''')
    
    Flow (prefix+'vp_noise',[prefix+'phi_noise',prefix+'shale'],
          '''
          math phi=${SOURCES[0]} shale=${SOURCES[1]}
          output="4090-960*shale-4510*phi"
          ''')
    
    Flow (prefix+'vs_noise',prefix+'vp_noise',
          'math vp=$SOURCE output="-880+0.782*vp"')

#
# Add correlated noise to the porosity in the reservoir interval
#

def add_noise_reservoir (prefix='',bk_std_dev=None,sd_std_dev=None):

    for i in range(15):
        Flow ('layer_%02d_noise' % (i+1),'layer_xypad sd_sim_%02d_real' % (i+1),
              'mask min=%g max=%g | dd type=float | add mode=m ${SOURCES[1]}' 
              % (i+0.5,i+1.5) )
    
    noises = []
    for i in range(15): noises.append('layer_%02d_noise' % (i+1))
        
    Flow ('sd_rfield',noises,'add mode=a ${SOURCES[1:15]}')
    
    sources = [prefix+'phi',prefix+'bk_sim_01_real',
               'ch_mask_xypad','sd_rfield',
               'as_mask_xypad','na_mask_xypad']
    
    Flow (prefix+'phi_noise',sources,
          '''
          math phi=${SOURCES[0]}
          bknoise=${SOURCES[1]} chmask=${SOURCES[2]}
          sdnoise=${SOURCES[3]} asmask=${SOURCES[4]} namask=${SOURCES[5]}
          output="phi+%g*bknoise*(1-chmask)+%g*sdnoise*(asmask+namask)"
          ''' % (bk_std_dev,sd_std_dev) )

#
# Add correlated noise to the porosity in the overburden interval
#

def add_noise_overburden (prefix='',bk_std_dev=None):

    Flow (prefix+'phi_noise',[prefix+'phi',prefix+'bk_sim_01_real'],
          '''
          math phi=${SOURCES[0]} bknoise=${SOURCES[1]}
          output="phi+%g*bknoise"
          ''' % (bk_std_dev) )

#
# Add padding in the x and y direction and reset the origin
#

def add_xy_pad (prefix='',xy_pad=0):

    Flow (prefix+'sand_xypad','sand',
          'pad beg1=%d beg2=%d end1=%d end2=%d | put o1=0 o2=0 o3=0'
          % (xy_pad,xy_pad,xy_pad,xy_pad) )

    for cube in ('layer','ch_mask','as_mask','na_mask'):
        Flow (cube+'_xypad',cube,
              'pad beg1=%d beg2=%d end1=%d end2=%d | put o1=0 o2=0 o3=0'
              % (xy_pad,xy_pad,xy_pad,xy_pad) )

#
# Apply a top taper to the property arrays
#

def apply_taper (   prefix='',memsize=100,
                    nz=None, dz=None, oz=None,
                    top_taper=None,     bot_taper=None,
                    top_taper_phi=None, top_taper_rho=None,
                    top_taper_vp=None,  top_taper_vs=None,
                    bot_taper_phi=None, bot_taper_rho=None,
                    bot_taper_vp=None,  bot_taper_vs=None):

    if (top_taper > 0):
        Flow (prefix+'top_taper',prefix+'bk_sim_01_real',
              'math output="(%g-x3)/%g" | clip2 upper=1'
              % ((nz-1)*dz,top_taper) )
    else:
        Flow (prefix+'top_taper',prefix+'bk_sim_01_real',
              'math output="1"')
        
    if (bot_taper > 0):
        Flow (prefix+'bot_taper',prefix+'bk_sim_01_real',
              'math output="x3/%g" | clip2 upper=1'
              % (bot_taper) )
    else:
        Flow (prefix+'bot_taper',prefix+'bk_sim_01_real',
              'math output="1"')
        
    taper_list = [prefix+'phi_noise',
                  prefix+'rho_noise',
                  prefix+'vp_noise',
                  prefix+'vs_noise']
        
    top_props = {prefix+'phi_noise'  : top_taper_phi,
                 prefix+'rho_noise'  : top_taper_rho,
                 prefix+'vp_noise'   : top_taper_vp,
                 prefix+'vs_noise'   : top_taper_vs}
    
    bot_props = {prefix+'phi_noise'  : bot_taper_phi,
                 prefix+'rho_noise'  : bot_taper_rho,
                 prefix+'vp_noise'   : bot_taper_vp,
                 prefix+'vs_noise'   : bot_taper_vs}
    
    for prop in (taper_list):
        Flow (prop+'_taper',[prop,prefix+'top_taper',prefix+'bot_taper'],
              '''
              math prop=${SOURCES[0]} top=${SOURCES[1]} bot=${SOURCES[2]}
              output="prop*top*bot+%g*(1-top)+%g*(1-bot)"
              ''' % (top_props[prop],bot_props[prop]) )

#
# Make reservoir model grids
#

def make_reservoir (    private=None, memsize=100, xy_pad=None,
                        nx=None,    ny=None,    nz=None,
                        dx=None,    dy=None,    dz=None,
                        ox=None,    oy=None,    oz=None,
                        oriu=None,  oriv=None,  oriw=None,
                        bk_ru=None, bk_rv=None, bk_rw=None,
                        sd_ru=None, sd_rv=None, sd_rw=None,
                        bk_std_dev=None,    sd_std_dev=None,
                        bd_depth=None,      md_depth=None,
                        bd_zshift=None,     md_zshift=None,
                        as_width=None,      as_shift=None, 
                        as_height0=None,    as_height1=None, 
                        as_shape=None,      as_sand=None,
                        bd_sand=None,       md_sand=None,
                        na_sand0=None,      na_sand1=None,
                        top_taper=None,     bot_taper=None,
                        top_taper_phi=None, top_taper_rho=None,
                        top_taper_vp=None,  top_taper_vs=None,
                        bot_taper_phi=None, bot_taper_rho=None,
                        bot_taper_vp=None,  bot_taper_vs=None):

# Get channel properties from the data server

    get_props(private=private)

# Put the channel properties into layer arrays

    make_prop_arrays(   memsize=memsize,
                        nx=nx, ny=ny,
                        dx=dx, dy=dy,
                        ox=ox, oy=oy)

# Make channel profile surfaces

    make_surfaces(  bd_depth=bd_depth,   md_depth=md_depth,
                    bd_zshift=bd_zshift, md_zshift=md_zshift)

# Erode the channel profile surfaces

    erode_surfaces (memsize=memsize)

# Put channel properties into a regular grid

    make_prop_grids (memsize=memsize,nz=nz,dz=dz,oz=oz)

# Make indicators and masks

    make_masks (memsize=memsize,nz=nz,dz=dz,oz=oz)

# Make sand fraction grid

    make_sand ( memsize=100,
                nz=nz, dz=dz, oz=oz,
                as_width=as_width,     as_shift=as_shift, 
                as_height0=as_height0, as_height1=as_height1, 
                as_shape=as_shape,     as_sand=as_sand,
                bd_sand=bd_sand,       md_sand=md_sand,
                na_sand0=na_sand0,     na_sand1=na_sand1)

# Add padding in the x and y direction and reset the origin

    add_xy_pad (prefix='res_',xy_pad=xy_pad)

# Make shale and porosity

    make_shale_and_phi (prefix='res_')
     
# Make correlated random fields 

    rfield.rfield(  name='res_bk_', seed=1, nr=1,
                    nx=nx+2*xy_pad, ny=ny+2*xy_pad, nz=nz,
                    dx=dx,          dy=dy,          dz=dz,
                    oriu=oriu,      oriv=oriv,      oriw=oriw,
                    ru=bk_ru,       rv=bk_rv,       rw=bk_rw)
    
    rfield.rfield(  name='sd_', seed=2, nr=15,
                    nx=nx+2*xy_pad, ny=ny+2*xy_pad, nz=nz,
                    dx=dx,          dy=dy,          dz=dz,
                    oriu=oriu,      oriv=oriv,      oriw=oriw,
                    ru=sd_ru,       rv=sd_rv,       rw=sd_rw)

# Add correlated noise to the porosity

    add_noise_reservoir (prefix='res_',
                         bk_std_dev=bk_std_dev,sd_std_dev=sd_std_dev)

# Make density and sonic velocities

    make_rho_and_vel (prefix='res_')

# Apply a top taper to the property arrays

    apply_taper (   prefix='res_',memsize=memsize,
                    nz=nz, dz=dz, oz=oz,
                    top_taper=top_taper,            bot_taper=bot_taper,
                    top_taper_phi=top_taper_phi,    top_taper_rho=top_taper_rho,
                    top_taper_vp=top_taper_vp,      top_taper_vs=top_taper_vs,
                    bot_taper_phi=bot_taper_phi,    bot_taper_rho=bot_taper_rho,
                    bot_taper_vp=bot_taper_vp,      bot_taper_vs=bot_taper_vs)

#
# Make overburden model grids
#

def make_overburden (   memsize=100, xy_pad=None,
                        nx=None,    ny=None,    nz=None,
                        dx=None,    dy=None,    dz=None,
                        ox=None,    oy=None,    oz=None,
                        oriu=None,  oriv=None,  oriw=None,
                        bk_ru=None, bk_rv=None, bk_rw=None,
                        bk_std_dev=None,
                        top_taper=None,     bot_taper=None,
                        top_taper_phi=None, top_taper_rho=None,
                        top_taper_vp=None,  top_taper_vs=None,
                        bot_taper_phi=None, bot_taper_rho=None,
                        bot_taper_vp=None,  bot_taper_vs=None):


# Make sand fraction grid

    grid = '''
       n1=%d d1=%g o1=0
       n2=%d d2=%g o2=0
       n3=%d d3=%g o3=0
       ''' % (nx+2*xy_pad,dx,ny+2*xy_pad,dy,nz,dz)

    Flow ('ovr_sand_xypad','','math'+grid+'output="0"',stdin=0)

# Make shale and porosity

    make_shale_and_phi (prefix='ovr_')
     
# Make correlated random fields 

    rfield.rfield(  name='ovr_bk_', seed=1, nr=1,
                    nx=nx+2*xy_pad, ny=ny+2*xy_pad, nz=nz,
                    dx=dx,          dy=dy,          dz=dz,
                    oriu=oriu,      oriv=oriv,      oriw=oriw,
                    ru=bk_ru,       rv=bk_rv,       rw=bk_rw)
    
# Add correlated noise to the porosity

    add_noise_overburden (prefix='ovr_',bk_std_dev=bk_std_dev)

# Make density and sonic velocities

    make_rho_and_vel (prefix='ovr_')

# Apply a top taper to the property arrays

    apply_taper (   prefix='ovr_',memsize=memsize,
                    nz=nz, dz=dz, oz=oz,
                    top_taper=top_taper,            bot_taper=bot_taper,
                    top_taper_phi=top_taper_phi,    top_taper_rho=top_taper_rho,
                    top_taper_vp=top_taper_vp,      top_taper_vs=top_taper_vs,
                    bot_taper_phi=bot_taper_phi,    bot_taper_rho=bot_taper_rho,
                    bot_taper_vp=bot_taper_vp,      bot_taper_vs=bot_taper_vs)
