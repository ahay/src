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

from rsf.proj import *
from math    import pi
import rfield
import os

# 
# Input parameter explanations.
#

#
# Most of the input parameters for make_reservoir and make_overburden
# are loaded into dictionaries:
#
# def make_reservoir (  memsize=100,            transpose memory size
#                       private=None,           data server parameters
#                       grid_par=None,          grid parameters
#                       geo_par=None,           geometrical parameters
#                       sand_par=None,          sand fraction parameters
#                       bk_noise_par=None,      background noise parameters
#                       sd_noise_par=None,      sand noise parameters
#                       taper_par=None):        taper parameters
#
# def make_overburden ( memsize=100,            transpose memory size
#                       grid_par=None,          grid parameters
#                       noise_par=None,         background noise parameters
#                       taper_par=None):        taper parameters
#
# Grid parameters
#
# The z direction parameters may be different for the reservoir and
# overburden, but for seismic modeling the x and y direction parameters
# should be the same.
#
# grid_par = {  'nx', 'dx', 'ox',   x grid parameters
#               'ny', 'dy', 'oy',   y grid parameters
#               'nz', 'dz', 'oz',   z grid parameters
#               'xy_pad'            x & y direction padding
#            }
#
# Notes on specification:
#
# This current version of channels.py interpolates 15 2D arrays of
# channel position and skew from Mathematica-exported files to your
# desired grid.  The original files were generated in Mathematica on
# the following grid:
#
# n1=181   d1=30m   o1=0m 
# n2=101   d2=30m   o2=0m          
#
# The Mathematica grid spans:
#
# 0m <= x <= 5400m  and  0m <= y <= 3000m
#
# Your grid should be within this range to avoid extrapolation on the
# edges.  These grid parameters set the location and size of your grid
# relative to the Mathematica arrays.  After the channels are created
# optional x and y padding is added, increasing nx and ny to values
# larger than those specified here (see padding parameter below), and
# the grid origin is reset to zero.
#
# Future versions of channels.py may generate position and skew
# internally, avoiding the interpolation-extrapolation issue.
#
# The following example spans:
#
# 0m <= x <= (nx-1)*dx = 5385m  and  0m <= y <= (ny-1)*dy = 2985m
#
# nx=360   dx=15m   ox=0m         
# ny=200   dy=15m   oy=0m          
#
# The z coordinate measures height (not depth) above the bottom of the
# first channel.  The top of the last channel is at 110m.
#
# There are two sets of vertical grid parameters, one for a reservoir
# grid that contains the channels, and another set for a grid of
# overburden with only background material with it's noise, but no
# channels.  The two grids should probably have the same x and y grid
# parameters, but they may have different z direction grid size,
# spacing, and location.
#
# You can choose any vertical grid.  Choosing a reservoir grid smaller
# than the channel range will produce a grid that clips some of the
# channels. Choosing a larger grid will generate extra "pad" volume
# above or below the channels.  After the channels are created the grid
# origin is reset to zero.  The overburden grid will have no channels,
# even if you choose it's location within the zone where the channels
# exist.
#
# The following example creates a vertical reservoir grid that spans
# -20m <= z <= 129m, with a 20m pad below the bottom channel, and a 19m
# pad above the top channel:
#
# nz=150 dz=1m oz=-20m
#
# Additional padding in the x and y directions: xy_pad
#
# This parameter controls the width (in cells) for additional all-shale
# padding to be added to all four horizontal sides of the model after
# the channels are constructed.  After the pad is added, the x and y
# origin is reset to zero.
#
# Channel geometrical properties
#
# geo_par = {                   amalgamated sand geometrical properties
#               'as_width'      width at channel bottom
#               'as_shift'      lateral shift at channel max bend
#               'as_height0'    height at channel inflection
#               'as_height1'    height at channel max bend
#               'as_shape'      cross-section shape factor
#
#                               bar drape geometrical properties
#               'bd_depth'      bar drape profile depth
#               'bd_zshift'     bar drape profile z shift
#
#                               margin drape geometrical properties
#               'md_depth'      margin drape profile depth
#               'md_zshift'     margin drape profile z shift
#           }
# 
# Sand fraction parameters 
#
# sand_par = {  'bd_sand'   bypass drape sand fraction (sf)
#               'md_sand'   margin drape sf
#               'as_sand'   amalgamated sand sf
#               'na_sand0'  non-amalgamated sand sf at channel top
#               'na_sand1'  non-amalgamated sand sf at amalgamated sand contact
#            }
#
# Porosity noise parameters
#
# The reservoir requires two of noise parameters, bk_noise_par and
# sd_noise_par, for the background and the sand.  A third set of noise
# parameters may be used for the overburden, although it commonly is
# the same as bk_noise_par.  See rfield.py for notes on the covariance
# parameters for random field generation.
#
# These parameters are used to generate random fields of noise to be
# added to the background and sand porosity.  The same noise propogates
# through the subsequent generation of the rho, vp, and vs grids.
#
# The code will make one realizataion for the background noise, and 15
# different realizations for the channels, one for each of the 15
# channels in the model.  The 15 sand realizatins are all made with the
# same set of parameters.  The random number seeds are preset.  A small
# coding change would allow for user setting of the random number seeds
# for different realizations of the entire model.
#
# noise_par = { 'taper_switch'  covariance taper switch
#               'std_dev'       porosity noise standard devietion
#               'alpha',        covariance shape parameter
#               'oriu',         covariance range orientation vectors
#               'oriv', 
#               'oriw',
#               'ru',           covariance range parameters
#               'rv',   
#               'rw'}
#
# Taper parameters
#
# A linear taper is applied, if desired, to the rock properties (phi,
# rho, vp, and vs) at the top and bottom of each of the reservoir and
# overburden grids.  However, it probably makes most sense to use only
# a top taper in the overburden, and only a bottom taper in the
# reservoir.
# 
# The properties are linearly interpolated between those generated by
# the channel model, and a set of constants for the top (or bottom) of
# the grid. The weighting is computed to produce the constant values at
# the top (or bottom) of the taper zone, the channel model output at
# the bottom (or top) of the taper zone and below (or above), and
# linear weighting with depth in between.
#
# If a taper thickness is set to zero no taper will be applied in the
# corresponding taper zone.
#
# taper_par = {             top taper parameters
#               'top_h'     thickness (m)
#               'top_phi'   porosity (fraction)
#               'top_rho'   density (gm/cc)
#               'top_vp'    Vp (m/s)
#               'top_vs'    Vs (m/s)
#
#                           bottom taper parameters
#               'bot_h'     thickness (m)
#               'bot_phi'   porosity (fraction)
#               'bot_rho'   density (gm/cc)
#               'bot_vp'    Vp (m/s)
#               'bot_vs'    Vs (m/s)
#             }

#
# Function definitions
#

#
# Convert props.dat format to RSF ascii
#

def format_props (target=None, source=None, env=None):

    inp   = open(str(source[0]),'r')    # open files
    file  = str(target[0])
    out   = open(file,'w')
    
    lines = inp.readlines()             # read input file
    lines.pop(0)                        # remove the first line
    for line in lines:                  # get all fields except the first one
        out.write(' '.join(line.split()[1:])+'\n')
        
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
        if not os.path.exists(type+'_data.rsf'):
            Fetch (type+'_data.rsf','jim',private)
    
    if not os.path.exists('props.dat'): # get props table from server
        Fetch ('props.dat','jim',private)
        
                                        # make a builder
    Props = Builder(action=Action(format_props),src_suffix='.dat',suffix='.asc')
    env   = Environment(BUILDERS={'Props':Props})
    
    env.Props('props')                  # run the builder(?)
                                        # convert ascii to native
    Flow ('props','props.asc','dd form=native')

#
# Put the channel properties into layer arrays
#

def make_prop_arrays (memsize=100, grid_par=None):
                        
    par = grid_par.copy()
    par['memsize'] = memsize

    for type in ('skew','pos'):         # interpolate skew and pos arrays
        Flow (type,type+'_data',
              '''
              dd form=native |
              remap1 n1=%(nx)d d1=%(dx)g o1=%(ox)g | 
              transp memsize=%(memsize)d plane=12  |
              remap1 n1=%(ny)d d1=%(dy)g o1=%(oy)g | 
              transp memsize=%(memsize)d plane=12
              ''' % (par) )

                                        # get the props out of the table
    Flow ('depth_list'  ,'props','window f1=0 n1=1')
    Flow ('wtop_list'   ,'props','window f1=1 n1=1')
    Flow ('wbot_list'   ,'props','window f1=2 n1=1')
    Flow ('aggrade_list','props','window f1=3 n1=1 | add add=20')

                                        # put the props in layer arrays
    for prop in ('depth','wtop','wbot','aggrade'):
        Flow (prop,prop+'_list',
              '''
              spray axis=2 n=%(ny)d d=%(dy)g o=%(oy)g | 
              spray axis=3 n=%(nx)d d=%(dx)g o=%(ox)g |
              transp memsize=%(memsize)d plane=13
              ''' % (par) )

#
# Make channel profile surfaces
#

def make_surfaces ( geo_par=None):

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
         
    prefix = ["ch_" ,"bd_"                 ,"md_"]
    depthf = [ 1    , geo_par['bd_depth']  , geo_par['md_depth']]
    zshift = [ 0    , geo_par['bd_zshift'] , geo_par['md_zshift']]
    
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
          rm -f bd_surface_15.rsf bd_erode_temp.rsf
          ''', stdout=0)

    Flow ('md_erode','md_surfaces ch_erode_temp',
          '''
          window squeeze=n f3=14 n3=1 > md_surface_15.rsf &&
          cat axis=3 ${SOURCES[1]} md_surface_15.rsf > md_erode_temp.rsf &&
          minmax mode=min file1=md_erode_temp.rsf file2=${SOURCES[0]} > $TARGET &&
          rm -f md_surface_15.rsf md_erode_temp.rsf
          ''', stdout=0)

#
# Put channel properties into a regular grid
#

def make_prop_grids (memsize=100, grid_par=None):

    par = grid_par.copy()
    par['memsize'] = memsize

    for prop in ('pos','skew','aggrade'):
        Flow (prop+'_pad',prop,'pad beg3=1 o3=0')
        Flow (prop+'_grid','ch_erode '+prop+'_pad',
              '''
              unif3 n1=%(nz)d d1=%(dz)g o1=%(oz)g layers=${SOURCES[1]} | 
              transp memsize=%(memsize)d plane=13 | 
              transp memsize=%(memsize)d plane=12
              ''' % (par) )
    
    for prop in ('depth','wtop','bd_erode','md_erode'):
        Flow (prop+'_pad' ,prop,
              'window squeeze=n f3=0 n3=1 | cat axis=3 $SOURCE')
        Flow (prop+'_grid','ch_erode '+prop+'_pad',
              '''
              unif3 n1=%(nz)d d1=%(dz)g o1=%(oz)g layers=${SOURCES[1]} | 
              transp memsize=%(memsize)d plane=13 | 
              transp memsize=%(memsize)d plane=12
              ''' % (par) )

#
# Make indicators and masks
#

def make_masks (memsize=100, grid_par=None):

    par = grid_par.copy()
    par['memsize'] = memsize

# layer indicator

    Flow ('layer','ch_erode',
          '''
          unif3 n1=%(nz)d d1=%(dz)g o1=%(oz)g 
          v00=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 |
          transp memsize=%(memsize)d plane=13 | 
          transp memsize=%(memsize)d plane=12
          ''' % (par) )
    
# above-layer-zero mask

    Flow ('layer_00_mask','layer','clip2 upper=1')
    
# channel mask

    Flow ('ch_mask','z_rel','mask max=-0.01 | dd type=float')

# bypass drape and margin drape masks

    if (par['nz'] > 1):
        Flow ('bd_mask','bd_erode_grid ch_mask',
              '''
              math surf=${SOURCES[0]} output="x3-surf" | 
              mask max=0 | dd type=float | add mode=m ${SOURCES[1]}
              ''')
        
        Flow ('md_mask','md_erode_grid ch_mask bd_mask',
              '''
              math surf=${SOURCES[0]} output="x3-surf" | 
              mask max=0 | dd type=float | add mode=m ${SOURCES[1]} |
              math bd=${SOURCES[2]} output="input*(input-bd)"
              ''')
    else:
        Flow ('bd_mask','bd_erode_grid ch_mask',
              '''
              math surf=${SOURCES[0]} output="%(oz)g-surf" | 
              mask max=0 | dd type=float | add mode=m ${SOURCES[1]}
              ''' % (par) )
              
        Flow ('md_mask','md_erode_grid ch_mask bd_mask',
              '''
              math surf=${SOURCES[0]} output="%(oz)g-surf" | 
              mask max=0 | dd type=float | add mode=m ${SOURCES[1]}
              math bd=${SOURCES[2]} output="input*(input-bd)"
              ''' % (par) )

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

def make_sand (grid_par=None, geo_par=None, sand_par=None):

    geo_par_copy = geo_par.copy()
    geo_par_copy['as_heightd'] = geo_par['as_height1']-geo_par['as_height0']

    sand_par_copy = sand_par.copy()
    sand_par_copy['na_sandd'] = sand_par['na_sand1']-sand_par['na_sand0']

    if (grid_par['nz'] > 1):
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
              output="m*(%(oz)g-a)/d" |
              clip2 upper=0
              ''' % (grid_par) )
         
    Flow ('pos_shift','pos_grid wtop_grid skew_grid',
          '''
          math p=${SOURCES[0]} w=${SOURCES[1]} s=${SOURCES[2]}
          output="p/w+%(as_shift)g*s/2"
          ''' % (geo_par_copy) )
         
    Flow ('as_height','skew_grid',
          'math s=$SOURCE output="%(as_height0)g+%(as_heightd)g*s"'
          % (geo_par_copy) )
         
    Flow ('sand_parm','pos_shift z_rel as_height',
          '''
          math x=${SOURCES[0]} z=${SOURCES[1]} h=${SOURCES[2]}
          output="-z/(1-h+4^%(as_shape)g*h*(x^2/%(as_width)g^2)^%(as_shape)g)"
          ''' % (geo_par_copy) )
         
    Flow ('sand_parm_clip','sand_parm','clip2 upper=1')
    
    formula = '''
               bd*%(bd_sand)g+md*%(md_sand)g+as*%(as_sand)g
              +na*(%(na_sand0)g+%(na_sandd)g*p)
              ''' % (sand_par_copy)
    
    Flow ('sand','sand_parm_clip bd_mask md_mask as_mask na_mask',
          '''
          math p=${SOURCES[0]} 
          bd=${SOURCES[1]} md=${SOURCES[2]} as=${SOURCES[3]} na=${SOURCES[4]}
          output="%s"
          ''' % (formula) )
     
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

def add_noise_reservoir (prefix='',bk_std_dev=0,sd_std_dev=0):

    for i in range(15):
        Flow ('layer_%02d_noise' % (i+1),'layer_xypad sd_sim_%02d_real' % (i+1),
              'mask min=%g max=%g | dd type=float | add mode=m ${SOURCES[1]}' 
              % (i+0.5,i+1.5) )
    
    noises = []
    for i in range(15): noises.append('layer_%02d_noise' % (i+1))
    
    Flow ('sd_rfield',noises,
          'cat axis=4 ${SOURCES[1:15]} | stack axis=4')
    
    sources = [prefix+'phi',prefix+'bk_sim_01_real',
               'ch_mask_xypad','sd_rfield',
               'as_mask_xypad','na_mask_xypad']
    
    Flow (prefix+'phi_noise',sources,
          '''
          math phi=${SOURCES[0]}
          bknoise=${SOURCES[1]} chmask=${SOURCES[2]}
          sdnoise=${SOURCES[3]} asmask=${SOURCES[4]} namask=${SOURCES[5]}
          output="phi+%g*bknoise*(1-chmask)"
          ''' % (bk_std_dev) )

#
# Add correlated noise to the porosity in the overburden interval
#

def add_noise_overburden (prefix='',bk_std_dev=0):

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
# Apply top & bottom tapers to the property arrays
#

def apply_taper (prefix='', memsize=100, grid_par=None, taper_par=None):
                    
#    grid_par = {'nz':nz, 'dz':dz}

    if (taper_par['top_h'] > 0):
        Flow (prefix+'top_taper',prefix+'bk_sim_01_real',
              'math output="(%g-x3)/%g" | clip2 upper=1'
              % ((grid_par['nz']-1)*grid_par['dz'],taper_par['top_h']) )
    else:
        Flow (prefix+'top_taper',prefix+'bk_sim_01_real',
              'math output="1"')
        
    if (taper_par['bot_h'] > 0):
        Flow (prefix+'bot_taper',prefix+'bk_sim_01_real',
              'math output="x3/%g" | clip2 upper=1'
              % (taper_par['bot_h']) )
    else:
        Flow (prefix+'bot_taper',prefix+'bk_sim_01_real',
              'math output="1"')
        
    taper_list = [prefix+'phi_noise',
                  prefix+'rho_noise',
                  prefix+'vp_noise',
                  prefix+'vs_noise']
        
    top_props = {prefix+'phi_noise'  : taper_par['top_phi'],
                 prefix+'rho_noise'  : taper_par['top_rho'],
                 prefix+'vp_noise'   : taper_par['top_vp'],
                 prefix+'vs_noise'   : taper_par['top_vs']}
    
    bot_props = {prefix+'phi_noise'  : taper_par['bot_phi'],
                 prefix+'rho_noise'  : taper_par['bot_rho'],
                 prefix+'vp_noise'   : taper_par['bot_vp'],
                 prefix+'vs_noise'   : taper_par['bot_vs']}
    
    for prop in (taper_list):
        Flow (prop+'_taper',[prop,prefix+'top_taper',prefix+'bot_taper'],
              '''
              math prop=${SOURCES[0]} top=${SOURCES[1]} bot=${SOURCES[2]}
              output="prop*top*bot+%g*(1-top)+%g*(1-bot)"
              ''' % (top_props[prop],bot_props[prop]) )

#
# Make reservoir model grids
#

def make_reservoir (    memsize=100, 
                        private=None, 
                        grid_par=None,
                        geo_par=None,
                        sand_par=None,
                        bk_noise_par=None,
                        sd_noise_par=None,
                        taper_par=None):

# Get channel properties from the data server

    get_props (private)

# Put the channel properties into layer arrays

    make_prop_arrays (memsize,grid_par)

# Make channel profile surfaces

    make_surfaces (geo_par)

# Erode the channel profile surfaces

    erode_surfaces (memsize)

# Put channel properties into a regular grid

    make_prop_grids (memsize,grid_par)

# Make indicators and masks

    make_masks (memsize,grid_par)

# Make sand fraction grid

    make_sand (grid_par,geo_par,sand_par)

# Add padding in the x and y direction and reset the origin

    add_xy_pad ('res_',grid_par['xy_pad'])

# Make shale and porosity

    make_shale_and_phi ('res_')
     
# Make correlated random fields 

    noise_grid_par = grid_par.copy()
    noise_grid_par['nx'] = grid_par['nx']+2*grid_par['xy_pad']
    noise_grid_par['ny'] = grid_par['ny']+2*grid_par['xy_pad']
                      
    real_par_bk = {'name':'res_bk_', 'nr':1,  'seed':1}
    real_par_sd = {'name':'sd_',     'nr':15, 'seed':2}
    
    rfield.rfield (real_par_bk,noise_grid_par,bk_noise_par)
    rfield.rfield (real_par_sd,noise_grid_par,sd_noise_par)

# Add correlated noise to the porosity

    add_noise_reservoir ('res_',bk_noise_par['std_dev'],sd_noise_par['std_dev'])

# Make density and sonic velocities

    make_rho_and_vel ('res_')

# Apply top & bottom tapers to the property arrays

    apply_taper ('res_',memsize,grid_par,taper_par)
    
#
# Make overburden model grids
#

def make_overburden (memsize=100, grid_par=None, 
                     noise_par=None, taper_par=None):

    noise_grid_par = grid_par.copy()
    noise_grid_par['nx'] = grid_par['nx']+2*grid_par['xy_pad']
    noise_grid_par['ny'] = grid_par['ny']+2*grid_par['xy_pad']
                      
# Make sand fraction grid

    grid = '''
           n1=%(nx)d d1=%(dx)g o1=0
           n2=%(ny)d d2=%(dy)g o2=0
           n3=%(nz)d d3=%(dz)g o3=0
           ''' % (noise_grid_par)

    Flow ('ovr_sand_xypad','','math'+grid+'output="0"',stdin=0)

# Make shale and porosity

    make_shale_and_phi ('ovr_')
     
# Make correlated random fields 

    real_par = {'name':'ovr_bk_', 'nr':1, 'seed':1}
    
    rfield.rfield (real_par,noise_grid_par,noise_par)

# Add correlated noise to the porosity

    add_noise_overburden ('ovr_',noise_par['std_dev'])

# Make density and sonic velocities

    make_rho_and_vel ('ovr_')

# Apply top & bottom tapers to the property arrays

    apply_taper ('ovr_',memsize,noise_grid_par,taper_par)
