#
# Compute spatial statistics in 1, 2, or 3 dimensions.
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
# February 2006
#
# $Id$
#

#
# Import RSF library.
#

from rsf.proj import *

#
# Fast spatial statistics computation with FFTs:
#   semi-variogram
#   general relative semi-variogram
#   covariance
#   correlogram
#
# Definitions for these spatial statistics can be found in:
#
# Deutsch, C. V., and Journel, A. G., 1998, GSLIB-Geostatistical
#   software library and user's guide: New York, Oxford University Press,
#   pp 44-45.
#
# Algorithm inspired from:
#
# Marcotte, D., 1996, Fast variogram computation with FFT, 
#   Computers and Geosciences, v. 22, n. 10, pp. 1175-1186.
# 

# 
# Input parameter explanations.
#

#
# The input parameters are loaded into a dictionary 'par':
#
# par = {   'ind', 'val':       Indicator and data filenames
#           'nx', 'ny', 'nz':   Array sizes
#           'dx', 'dy', 'dz':   Cell sizes
#       }
#

#
# Value and indicator data files: par['ind','val']
#
# These functions will compute a spatial statistics grid from the data
# values stored in the file named in par['val'].  The data must be
# stored in a regular 1, 2, or 3 dimensional array with array and cell
# sizes specified elsewhere in the dictionary "par" (see below). 
# However, missing values are allowed.
#
# The array named in par['ind'] should contain only zeros and ones
# indicating the absence (zero) or presence (one) of data in the
# corresponding value array.  The value array should have zero values
# in locations where data is missing.  The value and indicator arrays
# should share the same dimensions.
#

#
# Array sizes: par['nx','ny','nz']
#
# The code works for 1D, 2D, or 3D data grids.  For 2D or 1D set one
# or two of the array dimensions to 1.
#
# Use only even values (or 1) for nx, ny, and nz.  I'm not convinced
# the code is working correctly for odd values.
#
# The FFTs will be most efficient if the array sizes can be factored
# into combinations of 2, 3, 4, and 5.  Nevertheless, the FFTs may
# still be fast if a remaining factor bigger than 5 is not too large.
#
# Use integers.
#

#
# Cell sizes: par['dx','dy','dz']
#
# You may use integers or floats, conversion to float is done for
# you when necessary.
#

# 
# End of input parameter explanations.
#

# 
# Array parameter computations, rules, and target names.
#

def paramaters (par,rules):

#
# Compute padded (double), augmented (plus one), and half grid dimensions 
#

    if (par['nx'] > 1):
        par['nx_pad']  = par['nx']*2
        par['nx_aug']  = par['nx']+1
        par['nx_half'] = par['nx']/2
        par['nx_flip'] = 2
    else:
        par['nx_pad']  = par['nx']
        par['nx_aug']  = par['nx']
        par['nx_half'] = par['nx']
        par['nx_flip'] = 1

    if (par['ny'] > 1):
        par['ny_pad']  = par['ny']*2
        par['ny_aug']  = par['ny']+1
        par['ny_half'] = par['ny']/2
        par['ny_flip'] = 2
    else:
        par['ny_pad']  = par['ny']
        par['ny_aug']  = par['ny']
        par['ny_half'] = par['ny']
        par['ny_flip'] = 1

    if (par['nz'] > 1):
        par['nz_pad']  = par['nz']*2
        par['nz_aug']  = par['nz']+1
        par['nz_half'] = par['nz']/2
        par['nz_flip'] = 2
    else:
        par['nz_pad']  = par['nz']
        par['nz_aug']  = par['nz']
        par['nz_half'] = par['nz']
        par['nz_flip'] = 1

#
# Compute grid origin 
#

    par['ox'] = -par['nx_half']*par['dx']
    par['oy'] = -par['ny_half']*par['dy']
    par['oz'] = -par['nz_half']*par['dz']

#
# Set up forward and inverse FFT rules 
#

    fft_rule =      '''
                    pad n1=%(nx_pad)d n2=%(ny_pad)d n3=%(nz_pad)d | 
                    rtoc | 
                    fft3 axis=1 pad=1 |
                    fft3 axis=2 pad=1 |
                    fft3 axis=3 pad=1
                    ''' % par

    inv_fft_rule =  '''
                    fft3 axis=3 pad=1 inv=y | 
                    fft3 axis=2 pad=1 inv=y | 
                    fft3 axis=1 pad=1 inv=y | 
                    real
                    '''
                    
    rotate_rule =   'rotate rot1=%(nx)d rot2=%(ny)d rot3=%(nz)d ' % par
    
    window_rule =   '''
                    window f1=%(nx_half)d n1=%(nx_aug)d
                           f2=%(ny_half)d n2=%(ny_aug)d
                           f3=%(nz_half)d n3=%(nz_aug)d | 
                    put o1=%(ox)g o2=%(oy)g o3=%(oz)g
                        d1=%(dx)g d2=%(dy)g d3=%(dz)g
                    ''' % par
                    
    inv_rule = inv_fft_rule + '|' + rotate_rule + '|' + window_rule
#
# Set up common rules and target names
#

# Target names

    ind           = par['ind']
    val           = par['val']

    ind_fft       = ind + '_fft'
    val_fft       = val + '_fft'
    val_sqr       = val + '_sqr'
    val_sqr_fft   = val_sqr + '_fft'

    par['pairs']  = ind + '_pairs'
    par['sumvar'] = val + '_sumvar'
    par['sumvv']  = val + '_sumvv'
    par['sumiv']  = val + '_sumiv'
    par['sumiv2'] = val + '_sumiv2'
    par['sumvi']  = val + '_sumvi'
    par['sumv2i'] = val + '_sumv2i'

# Data squaring rule

    rules[val_sqr] = [val,'math output="input^2"']

# FFT rules

    rules[ind_fft]     = [ind,    fft_rule]
    rules[val_fft]     = [val,    fft_rule]
    rules[val_sqr_fft] = [val_sqr,fft_rule]

# Sum rules

    rules[par['pairs']]  = [ind_fft,
    'math output="input*conj(input)" | %s' % inv_rule]

    rules[par['sumvv']]  = [val_fft,
    'math output="input*conj(input)" | %s' % inv_rule]

    rules[par['sumiv']]  = [[ind_fft,val_fft],
    'math output="i*conj(v)" i=${SOURCES[0]} v=${SOURCES[1]} | %s' % inv_rule]

    rules[par['sumiv2']] = [[ind_fft,val_sqr_fft],
    'math output="i*conj(v)" i=${SOURCES[0]} v=${SOURCES[1]} | %s' % inv_rule]

    rules[par['sumvi']]  = [[ind_fft,val_fft],
    'math output="v*conj(i)" i=${SOURCES[0]} v=${SOURCES[1]} | %s' % inv_rule]

    rules[par['sumv2i']] = [[ind_fft,val_sqr_fft],
    'math output="v*conj(i)" i=${SOURCES[0]} v=${SOURCES[1]} | %s' % inv_rule]

    rules[par['sumvar']] = [[ind_fft,val_fft,val_sqr_fft],
        '''
        math output="(i*conj(v2)+v2*conj(i))/2-v*conj(v)"
        i=${SOURCES[0]} v=${SOURCES[1]} v2=${SOURCES[2]} | %s
        ''' % inv_rule]

#
# Variogram rules
#

def variogram (par,rules):

    paramaters(par,rules)

    rules[par['val']+'_var'] = [[par['pairs'],par['sumvar']],
        'clip2 lower=1 | math output="sumvar/input" sumvar=${SOURCES[1]}']

#
# General relative variogram rules
#

def gr_variogram (par,rules):

    paramaters(par,rules)

    sources = [par['pairs'],par['sumvar'],par['sumiv'],par['sumvi']]

    rules[par['val']+'_vgr'] = [sources,
        '''
        math output="4*input*sumvar/((sumiv+sumvi)^2)"
        sumvar=${SOURCES[1]} sumiv=${SOURCES[2]} sumvi=${SOURCES[3]}
        ''']

#
# Covariance rules
#

def covariance (par,rules):

    paramaters(par,rules)
    
    sources = [par['pairs'],par['sumvv'],par['sumiv'],par['sumvi']]

    rules[par['val']+'_cov'] = [sources,
        '''
        clip2 lower=1 | 
        math output="(sumvv-sumiv*sumvi/input)/input"
        sumvv=${SOURCES[1]} sumiv=${SOURCES[2]} sumvi=${SOURCES[3]}
        ''']

#
# Correlogram rules
#

def correlogram (par,rules):

    paramaters(par,rules)

    sources = [par['pairs'],par['sumvv'],
               par['sumiv'],par['sumiv2'],
               par['sumvi'],par['sumv2i']]

    corr_rule = """
                     (sumvv -sumiv*sumvi/input)/
                sqrt((sumiv2-sumiv*sumiv/input)*(sumv2i-sumvi*sumvi/input))
                """

    rules[par['val']+'_cor'] = [sources,
        '''
        clip2 lower=1 | 
        math output="%s"
        sumvv=${SOURCES[1]} 
        sumiv=${SOURCES[2]} sumiv2=${SOURCES[3]}
        sumvi=${SOURCES[4]} sumv2i=${SOURCES[5]}
        ''' % (corr_rule)]
