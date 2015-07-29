#
# Generate spatially correlated Gaussian random fields in 1, 2, or 3
# dimensions, with a possibly anisotropic, skewed, and rotated "stable"
# covariance model.  The generated Gaussian random fields will have a
# mean of zero and a variance of one.
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
# $Id: rfield.py 5879 2010-05-03 01:44:26Z ivlad $
#

#
# Import RSF and math libraries.
#

from rsf.proj import *
from math    import sqrt,pi

def rfield (real_par,grid_par,covar_par):

# 
# Input parameter explanations.
#

# The input parameters are loaded into dictionaries:
#
# real_par = {  'name':   File name prefix
#               'nr':     Number of realizations
#               'seed':   Random number seed
#            }
#
# grid_par = {  'nx', 'ny', 'nz':   Array sizes
#               'dx', 'dy', 'dz':   Cell sizes
#            }
#
# covar_par = {'taper_switch':      Covariance taper switch
#              'alpha':             Covariance shape parameter
#              'oriu':[ 5,1,0]      Covariance range orientation vectors
#              'oriv':[-1,5,0] 
#              'oriw':[ 0,0,1]
#              'ru'  :              Covariance range parameters
#              'rv'  :
#              'rw'  :
#             }
#

#
# File name prefix: real_par['name']
#
# A string to prepend to each file generated.
#

#
# Number of realizations to generate: real_par['nr']
#
# Use an integer in the range 1 <= nr <= 99.
#

#
# Random number seed: real_par['seed']
#
# The random number seed initializes the pseudo random number
# generator.  Repeat runs of the script with the same random number
# seed will produce exactly the same random field realization, except
# perhaps for roundoff differences, even on different platforms.  For
# different realizations choose different random number seeds.
#
# Use an integer.
#

#
# Array sizes: grid_par['nx','ny','nz']
#
# The code works for 1D, 2D, or 3D simulations.  For 2D or 1D set one
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
# Cell sizes: grid_par['dx','dy','dz']
#
# You may use integers or floats, conversion to float is done for
# you when necessary.
#

#
# Covariance grid cosine taper switch: covar_par['taper_switch']
#
# A cosine taper will be applied to the edges of the grid if
# taper_switch is true, otherwise not.
#
# The default is to use a cosine taper, but some covariance models work
# better without it.  An example of a covariance model that works
# better without a taper is the stable model with small alpha and with
# grid-aligned anisotropy.  When alpha is small the covariance can
# decay very slowly so that it still has significant values, but small
# derivatives on the edges of the grid. Because the derivatives are
# small a taper is not needed.  However, if the cosine taper is applied
# anyway, the covariance will end up with an unintended dominant cosine
# shape on the edges, producing unexpected results.
#
# Anisotropic covariance models not aligned with the grid are more
# difficult.  These models have a discontinuity on the grid edges
# because of the spatial periodicity of the FFT which is aligned with
# the grid, not the anisotropy.  This discontinuity will also produce
# grid-aligned linear artifacts.  I don't know a good way to make
# simulations with small alpha and non-grid-aligned anisotropy yet.
# There are some situations where it seems you must put up with one
# kind of artifact or another depending on whether tapering is turned
# on or not.  It might help to apply averaging at the discontinuities,
# but I haven't tried that yet. (see Briggs and Henson, The DFT: An
# Owner's Manual for the Discrete Fourier Transform, SIAM, look in the
# index for AVED, "average values at endpoints and discontinuities.")
#
# Some other time I might think of a better way of tapering that is
# more robust and requires less user knowledge...
#

#
# Covariance range orientation vectors: covar_par['oriu','oriv','oriw']
#
# The orientation vectors define a possibly rotated and skewed
# coordinate system orienting the covariance model.  The correlation
# ranges will be specified along these vectors.
#
# The orientation vectors do not need to be unit length, they will be
# scaled.  The orientation vectors do not need to be orthogonal,
# allowing for skewed covariance models.
#
# You may use integers or floats, conversion to float is done for you
# when necessary.
#
# Anisotropic, orthogonal, and rotated example
# oriu = [ 5,1,0]
# oriv = [-1,5,0]
# oriw = [ 0,0,1]
#
# Isotropic example
# oriu = [1,0,0]
# oriv = [0,1,0]
# oriw = [0,0,1]
#

#
# Covariance range parameters: covar_par['ru','rv','rw']
#
# Each covariance range parameter controls the correlation range in the
# direction of the corresponding orientation vector.
# 
# The range parameters should be small enough so that the derivative of
# the covariance decays to "small" values at the edges of the grid.
# This reflects the intuitive idea that it does not make sense to model
# correlations that go beyond the edge of the grid.  The edges of the
# covariance model can be tapered to prevent some of the artifacts that
# would otherwise be present, but range parameters that are too large
# will produce a covariance dominated by the taper rather than the
# desired model.
#
# You may specify integers or floats, conversion to float is done for
# you when necessary.
#

# 
# "Stable" covariance model shape parameter: covar_par['alpha']
#
# The stable covariance model for unit variance is:
#
# cov(h) = exp(-(h/r)^alpha)
#
# where h is a separation distance, or "lag," in a certain direction, r
# is the corresponding range parameter, and alpha is the shape
# parameter.
#
# The corresponding stable semi-variogram model is:
#
# gamma(h) = 1-exp(-(h/r)^alpha)
#
# The stable covariance model is asymptotically fractal (power law)
# for small enough h:
#
# cov(h->small) = (h/r)^alpha
#
# For larger values of h, cov(h) approaches zero and gamma(h)
# approaches a "sill" of one.
#
# The shape parameter alpha must be in the range 0 <= alpha <= 2.  The
# special cases of 0, 1, and 2 correspond to the familiar "Nugget",
# "Exponential", and "Gaussian" models of geostatistics.  These models
# produce different amounts of "roughness" in the correlated field. The
# Nugget model produces rough spatially uncorrelated random fields that
# are discontinuous everywhere (Gaussian white noise).  The Exponential
# model produces spatially correlated and continuous random fields with
# intermediate roughness and first derivatives that are discontinuous
# everywhere. The Gaussian model produces smooth spatially correlated
# random fields with continuous derivatives.  Fractional values of
# alpha are allowed, producing random fields anywhere in this spectrum
# of roughness.
# 
# You may specify an integer or a float, conversion to float is done
# for you when necessary.
#

# 
# End of input parameter explanations.
#

#
# Compute the grid origin.
#
# To construct the random field we will need to sample the covariance
# model on a grid and compute its FFT (among other things).  The FFT
# expects spatially periodic data with "wrap-around" storage.  That is,
# the origin of the coordinate system, (x,y,z)=(0,0,0), is expected to
# coincide with the corners of the grid, and the parts of the
# covariance model corresponding to negative values of x, y, or z are
# expected to be stored in the "second half" of the x, y, or z axis.
#
# However, this storage order is non-intuitive and inconvenient for
# sampling the covariance and taper functions.
#
# So instead we will set the values of (x,y,z) at the grid origin so
# that (0,0,0) is in the center of the array.  This results in more
# intuitive covariance and taper models centered in the middle of the
# grid and more convenient sampling.  Later, we will rotate the tapered
# covariance grid to shift the origin back into the corners of the
# grid. 
#

    grid_par['ox'] = -float(grid_par['nx']/2)*grid_par['dx']
    grid_par['oy'] = -float(grid_par['ny']/2)*grid_par['dy']
    grid_par['oz'] = -float(grid_par['nz']/2)*grid_par['dz']

#
# Set up the grid specification.
#

    grid_str = '''
               n1=%(nx)d d1=%(dx)g o1=%(ox)g
               n2=%(ny)d d2=%(dy)g o2=%(oy)g
               n3=%(nz)d d3=%(dz)g o3=%(oz)g
               label1="x" label2="y" label3="z"
               ''' % (grid_par)
               
    grid_par['string'] = grid_str
#
# Scale the orientation vectors.
#

    lenu = sqrt(covar_par['oriu'][0]**2
               +covar_par['oriu'][1]**2
               +covar_par['oriu'][2]**2)
               
    lenv = sqrt(covar_par['oriv'][0]**2
               +covar_par['oriv'][1]**2
               +covar_par['oriv'][2]**2)
               
    lenw = sqrt(covar_par['oriw'][0]**2
               +covar_par['oriw'][1]**2
               +covar_par['oriw'][2]**2)

    uniu = [float(covar_par['oriu'][0])/lenu,
            float(covar_par['oriu'][1])/lenu,
            float(covar_par['oriu'][2])/lenu]
            
    univ = [float(covar_par['oriv'][0])/lenv,
            float(covar_par['oriv'][1])/lenv,
            float(covar_par['oriv'][2])/lenv]
            
    uniw = [float(covar_par['oriw'][0])/lenw,
            float(covar_par['oriw'][1])/lenw,
            float(covar_par['oriw'][2])/lenw]

#
# Compute a scaled distance variable.
#
# This array measures the distance from the origin after scaling by a
# specified covariance range parameter along each corresponding
# orientation vector.  The scaling allows the user to control the
# covariance range and anisotropy in the possibly skewed and rotated
# uvw coordinate system.
#
# This calculation is easier to understand in two steps, transformation
# to the uvw coordinate system, and calculation of the scaled distance.
# Here they are combined into one step to avoid some unnecessary
# intermediate target files.
#
# The transformation step works like this:
#
# u = x*uniu[0] + y*uniu[1] + z*uniu[2]
# v = x*univ[0] + y*univ[1] + z*univ[2]
# w = x*uniw[0] + y*uniw[1] + z*uniw[2]
#
# The scaled distance calculation works like this:
#
# dist = sqrt((u/ru)^2 + (v/rv)^2 + (w/rw)^2)
#

    dist = real_par['name']+'dist'

    formula = '''
              output="sqrt(((x1*(%g)+x2*(%g)+x3*(%g))/(%g))^2
                          +((x1*(%g)+x2*(%g)+x3*(%g))/(%g))^2
                          +((x1*(%g)+x2*(%g)+x3*(%g))/(%g))^2)"
              ''' % (uniu[0],uniu[1],uniu[2],covar_par['ru'],
                     univ[0],univ[1],univ[2],covar_par['rv'],
                     uniw[0],uniw[1],uniw[2],covar_par['rw'])
                     
    Flow (dist,'','math'+grid_str+formula,stdin=0)
         
#
# Make a "stable" covariance grid.
#

    covar = real_par['name']+'covar'
    
    Flow (covar,dist,'math output="exp(-(input^%(alpha)g))"' % (covar_par))

#
# Make a cosine taper.
#
# The covariance should be smoothly tapered if necessary to eliminate
# any non-zero derivatives on the edges of the grid.  Otherwise the
# spatially periodic covariance will have a sharp "ridge" on the grid
# edge, which in turn produces a corresponding ridge on the weight
# function. Convolution of a random noise with a moving average weight
# function having a sharp ridge along the edges will produce linear
# artifacts throughout the grid aligned with the grid axes.
#
    
    taper = real_par['name']+'taper'
    
    formula = '''
              output="(cos(x1*%g)+1)*(cos(x2*%g)+1)*(cos(x3*%g)+1)/8"
              ''' % (2*pi/(grid_par['nx']*grid_par['dx']),
                     2*pi/(grid_par['ny']*grid_par['dy']),
                     2*pi/(grid_par['nz']*grid_par['dz']))

    if covar_par['taper_switch']:
        Flow (taper,'','math'+grid_str+formula,stdin=0)
    else:
        Flow (taper,'','spike'+grid_str,stdin=0)

#
# Apply the taper to the covariance, and rotate the coordinate system
# origin to the grid origin.
#

    rotate = real_par['name']+'rotate'

    rotx = int(grid_par['nx']/2)
    roty = int(grid_par['ny']/2)
    rotz = int(grid_par['nz']/2)
        
    Flow (rotate,[covar,taper],
          '''
          add mode=m ${SOURCES[1]} |
          rotate rot1=%d rot2=%d rot3=%d
          ''' % (rotx,roty,rotz) )

#
# Make uncorrelated Gaussian random noise with unit variance.
#

    noise = real_par['name']+'noise'
    
    Flow (noise,covar,
          '''
          put n4=%(nr)d d4=1 o4=1 label4="r" | 
          noise var=1 rep=y seed=%(seed)d
          ''' % (real_par) )
     
#
# Compute the FFT of the covariance and random noise.
#
# The Fourier transform of the covariance is the power spectrum. 
# Because the covariance should be real and symmetric (an even
# function), the power spectrum should have no imaginary component,
# except for roundoff.  If there is a significant imaginary component
# then the covariance does not have the proper symmetry, or has not
# been rotated correctly (anything else?).
#
# Also, the power spectrum should be everywhere positive, except for
# roundoff.  If not, the covariance is not a valid positive-definite
# model.  A non-positive-definite covariance implies a complex valued
# random field.
#

    pspec       = real_par['name']+'pspec'
    pspec_real  = real_par['name']+'pspec_real'
    noise_fft   = real_par['name']+'noise_fft'

    Flow (pspec,rotate,
          '''
          rtoc | 
          fft3 pad=1 opt=n axis=1 | 
          fft3 pad=1 opt=n axis=2 | 
          fft3 pad=1 opt=n axis=3 |
          put label1="wx" label2="wy" label3="wz"
          ''')
    
    Flow (pspec_real,pspec,'real')
    
    Flow (noise_fft,noise,
          '''
          rtoc | 
          fft3 pad=1 opt=n axis=1 | 
          fft3 pad=1 opt=n axis=2 | 
          fft3 pad=1 opt=n axis=3 |
          put label1="wx" label2="wy" label3="wz"
          ''')
    
# 
# Compute the amplitude spectrum.
# 
# The amplitude spectrum we will use here is simply the square root of
# the power spectrum.  This selects the unique real-valued amplitude
# spectrum from the family of possible complex-valued and
# Hermite-symmetric amplitude spectra.
#
# If the power spectrum has any negative values, corresponding to an
# invalid non-positive-definite covariance model, then the amplitude
# spectrum will end up complex valued, but not Hermite symmetric,
# implying a complex-valued random field.
# 
# The stable covariance model is positive definite for 0 <= alpha < 2.
# The Gaussian model, alpha=2 is on the boundary of positive
# definiteness.  It is not positive definite, but it is non-negative
# definite.  Values larger than 2 produce negative definite models.
# The use of a Gaussian convariance model with a discrete Fourier
# transform will produce some small nagative values due to round off,
# but these can be clipped without introducing significant artifacts.
#
# There are many possible ways we could deal with imaginary and/or
# negative values in the power spectrum.  Here we just discard the
# imaginary part and clip any real values below zero.  Thus, there
# is nothing to stop the user from specifing a non-positive-definite
# covariance, but if he does it will truncated in Fourier space
# into something that is positive definite.
#
# Additional experimentation might turn up better methods for dealing
# with complex or negative power spectrum values.
#

    aspec      = real_par['name']+'aspec'
    aspec_real = real_par['name']+'aspec_real'

    Flow (aspec_real,pspec_real,'clip2 lower=0 | add sqrt=1')
    Flow (aspec,     aspec_real,'rtoc')

#
# Combine the amplitude spectrum and noise.
#
# Multiplication in Fourier space is the same thing as convolution in
# real space.  Thus, this step is the same thing as performing a
# weighted moving average of the uncorrelated Gaussian noise, producing
# a correlated Gaussian noise.
#
# If we wanted to see the shape of the moving average weight function
# implied by this Fourier space convolution, we could calculate the
# inverse Fourier transform of the amplitude spectrum.  The weight
# function is sometimes called the "covariance square root" or
# "convolution square root" of the covariance function.  That is, the
# weight function is a function which when convolved with a reversed
# copy of itself produces the desired covariance function.
#
# Because we selected the real-valued amplitude spectrum, rather than a
# more general complex-valued and Hermite-symmetric amplitude spectrum,
# the weight function is symmetric.  Otherwise the weight function
# could be asymmetric, which might be interesting to tinker with some
# day...
#

    sim_i_fft = real_par['name']+'sim_%02d_fft'

    for i in range(real_par['nr']): 
        Flow (sim_i_fft % (i+1),[noise_fft,aspec],
              '''
              window squeeze=n f4=%d n4=1 | 
              add mode=m ${SOURCES[1]}
              ''' % (i) )

#
# Invert the simulation FFT.
#
# The resulting stochastic simulation will have both real and imaginary
# parts.  However, the imaginary part is produced by roundoff and should
# be small in magnitude.
#

    sim_i      = real_par['name']+'sim_%02d'
    sim_i_real = real_par['name']+'sim_%02d_real'

    for i in range(real_par['nr']): 
        Flow (sim_i % (i+1),sim_i_fft % (i+1),
              '''
              fft3 inv=y pad=1 opt=n axis=3 |
              fft3 inv=y pad=1 opt=n axis=2 |
              fft3 inv=y pad=1 opt=n axis=1 |
              put label1="x" label2="y" label3="z"
              ''')
    
        Flow (sim_i_real % (i+1),sim_i % (i+1),
              'real | put o1=0 o2=0 o3=0')
