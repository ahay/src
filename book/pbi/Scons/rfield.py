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
# $Id: rfield.py 594 2005-07-25 08:14:14Z jennings $
#

#
# Import RSF and math libraries.
#

from rsfproj import *
from math    import *

def rfield (name, 
            nr = 1,                 # number of realizations
            seed = 1,               # random number seed
            taper_switch = True,    # covariance-model taper switch
            nx = 1, ny = 1, nz = 1, # array sizes
            dx = 1, dy = 1, dz = 1, # cell sizes
            oriu = [1,0,0],         # covariance range orientation vectors
            oriv = [0,1,0],
            oriw = [0,0,1],
            ru = 1,                 # covariance range parameters
            rv = 1,
            rw = 1,
            alpha = 1):             # covariance shape parameter

# 
# Input parameter explanations.
#

#
# Array sizes: nx, ny, & nz.
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
# Number of realizations to generate: nr.
#
# Use an integer in the range 1 <= nr <= 99.
#

#
# Random number seed: seed.
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
# Covariance grid cosine taper switch: taper_switch.
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
# Cell sizes: dx, dy, dz.
#
# You may use integers or floats, conversion to float is done for
# you when necessary.
#

#
# Covariance range orientation vectors: oriu, oriv, oriw.
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
# Covariance range parameters: ru, rv, rw.
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
# "Stable" covariance model shape parameter: alpha.
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

    ox = -float(nx//2)*dx
    oy = -float(ny//2)*dy
    oz = -float(nz//2)*dz

#
# Set up the grid specification.
#

    grid = '''
       n1=%d d1=%g o1=%g
       n2=%d d2=%g o2=%g
       n3=%d d3=%g o3=%g
       ''' % (nx,dx,ox,ny,dy,oy,nz,dz,oz)

#
# Scale the orientation vectors.
#

    lenu = sqrt(oriu[0]**2+oriu[1]**2+oriu[2]**2)
    lenv = sqrt(oriv[0]**2+oriv[1]**2+oriv[2]**2)
    lenw = sqrt(oriw[0]**2+oriw[1]**2+oriw[2]**2)

    uniu = [float(oriu[0])/lenu,float(oriu[1])/lenu,float(oriu[2])/lenu]
    univ = [float(oriv[0])/lenv,float(oriv[1])/lenv,float(oriv[2])/lenv]
    uniw = [float(oriw[0])/lenw,float(oriw[1])/lenw,float(oriw[2])/lenw]

#
# Compute u, v, and w arrays.
#
# These arrays hold the u, v, and w coordinate values for a possibly
# skewed, and rotated coordinate system defined by the covariance range
# orientation vectors.
#

    u = name+'u'
    v = name+'v'
    w = name+'w'

    formula = grid+'output="x1*(%g)+x2*(%g)+x3*(%g)"'
    
    Flow (u,'','math'+formula % (uniu[0],uniu[1],uniu[2]),stdin=0)
    Flow (v,'','math'+formula % (univ[0],univ[1],univ[2]),stdin=0)
    Flow (w,'','math'+formula % (uniw[0],uniw[1],uniw[2]),stdin=0)

#
# Compute a scaled distance variable.
#
# This array measures the distance from the origin after scaling by a
# specified covariance range parameter along each corresponding
# orientation vector.  The scaling allows the user to control the
# covariance range and anisotropy in the possibly skewed and rotated
# uvw coordinate system.
#

    dist = name+'dist'

    Flow (dist,[u,v,w],
          '''
          math u=${SOURCES[0]} v=${SOURCES[1]} w=${SOURCES[2]}
          output="sqrt((u/%g)^2+(v/%g)^2+(w/%g)^2)"
          ''' % (ru,rv,rw),stdin=0
         )
         
#
# Make a "stable" covariance grid.
#

    covar = name+'covar'
    
    Flow (covar,dist,'math output="exp(-(input^%g))"' % (alpha))

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
    
    taper = name+'taper'
    
    if (taper_switch == 1):
        formula = grid+'output="(cos(x1*%g)+1)*(cos(x2*%g)+1)*(cos(x3*%g)+1)/8"'
        Flow (taper,'',
              'math'+formula % (2*pi/(nx*dx),2*pi/(ny*dy),2*pi/(nz*dz)),
              stdin=0
             )
    else:
        Flow (taper,'','spike'+grid,stdin=0)

#
# Apply the taper to the covariance.
#

    covar_taper = name+'covar_taper'

    Flow (covar_taper,[covar,taper],'add mode=m ${SOURCES[1]}')

#
# Shift the coordinate system origin to the grid origin.
#
# The shift is accomplished by cutting the grid in half and putting
# the halves back together in the reverse order, along each grid axis.
#

    shift_x = name+'shift_x'
    shift_y = name+'shift_y'
    shift   = name+'shift'

    if (nx > 1):
        Flow (shift_x,covar_taper,
              '''
              window n1=%d > half_x.rsf &&
              window f1=%d < $SOURCE | cat axis=1 half_x.rsf > $TARGET &&
              rm half_x.rsf
              ''' % (nx/2,nx/2), stdout=0
             )
    else:
        Flow (shift_x,covar_taper,
              'cp $SOURCE $TARGET',stdin=0,stdout=0
             )
    
    if (ny > 1):
        Flow (shift_y,shift_x,
              '''
              window n2=%d > half_y.rsf &&
              window f2=%d < $SOURCE | cat axis=2 half_y.rsf > $TARGET &&
              rm half_y.rsf
              ''' % (ny/2,ny/2), stdout=0
             )
    else:
        Flow (shift_y,shift_x,
              'cp $SOURCE $TARGET',stdin=0,stdout=0
             )
    
    if (nz > 1):
        Flow (shift,shift_y,
              '''
              window n3=%d > half_z.rsf &&
              window f3=%d < $SOURCE | cat axis=3 half_z.rsf > $TARGET &&
              rm half_z.rsf
              ''' % (nz/2,nz/2), stdout=0
             )
    else:
        Flow (shift,shift_y,
              'cp $SOURCE $TARGET',stdin=0,stdout=0
             )

#
# Make uncorrelated Gaussian random noise with unit variance.
#

    noise   = name+'noise'
    noise_i = name+'noise_%02d'
    
    Flow (noise,covar,
          'pad n4=%d | put n4=%d d4=1 o4=1 | noise var=1 rep=y seed=%d' 
          % (nr,nr,seed)
         )
     
    for i in range(nr):
        Flow (noise_i % (i+1),noise,'window f4=%d n4=1' % (i))
    
#
# Compute the FFT of the covariance and random noise.
#
# The Fourier transform of the covariance is the power spectrum. 
# Because the covariance should be real and symmetric (an even
# function), the power spectrum should have no imaginary component,
# except for roundoff.  If there is a significant imaginary component
# then the covariance does not have the proper symmetry, or has not
# been shifted correctly (anything else?).
#
# Also, the power spectrum should be everywhere positive, except for
# roundoff.  If not, the covariance is not a valid positive-definite
# model.  A non-positive-definite covariance implies a complex valued
# random field.
#

    pspec       = name+'pspec'
    pspec_real  = name+'pspec_real'
    noise_i_fft = name+'noise_%02d_fft'

    Flow (pspec,shift,
          'rtoc | fft3 pad=1 axis=1 | fft3 pad=1 axis=2 | fft3 pad=1 axis=3')
    
    Flow (pspec_real,pspec,'real')
    
    for i in range(nr): 
        Flow (noise_i_fft % (i+1),noise_i % (i+1),
              'rtoc | fft3 pad=1 axis=1 | fft3 pad=1 axis=2 | fft3 pad=1 axis=3')

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

    aspec      = name+'aspec'
    aspec_real = name+'aspec_real'

    Flow (aspec,pspec_real,'clip2 lower=0 | rtoc | add sqrt=1')
    Flow (aspec_real,aspec,'real')

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

    sim_i_fft = name+'sim_%02d_fft'

    for i in range(nr): 
        Flow (sim_i_fft % (i+1),[aspec,noise_i_fft % (i+1)],
              'add mode=m ${SOURCES[1]}'
             )

#
# Invert the simulation FFT.
#
# The resulting stochastic simulation will have both real and imaginary
# parts.  However, the imaginary part is produced by roundoff and should
# be small in magnitude.
#

    sim_i      = name+'sim_%02d'
    sim_i_real = name+'sim_%02d_real'

    for i in range(nr): 
        Flow (sim_i % (i+1),sim_i_fft % (i+1),
              '''
              fft3 inv=y pad=1 axis=3 |
              fft3 inv=y pad=1 axis=2 |
              fft3 inv=y pad=1 axis=1
              '''
             )
    
        Flow (sim_i_real % (i+1),sim_i % (i+1),'real')
