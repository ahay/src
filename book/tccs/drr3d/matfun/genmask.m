function mask = genmask(u, r, type,seed)
%GENMASK:Generate Random Sampling Mask
%
%	 mask = genmask(u,r,type)
%    u, image
%    r, data KNOWN ratio
%    type: data lose type
%   'r': random lose rows
%   'c': random lose columns
%   'p': random lose pixel
%   'seed': seed of random number generator
%
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details:
%  http://www.gnu.org/licenses/
%
%  References:   
%
%  [1] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  [2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  [3] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  [4] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  [5] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.

[m,n] = size(u);
mask = zeros(m,n);

switch type
    case 'r'
        row = rperm(m,seed);
        k = fix(r*m);
        row = row(1:k);
        mask(row,:) = 1;
    case 'c'
        column = rperm(n,seed);
        k = fix(r*n);
        column = column(1:k);
        mask(:, column) = 1;
    case 'p'
        pix = rperm(m*n,seed);
        r = fix(r*m*n);
        pix = pix(1:r);
        mask(pix) = 1;
end

