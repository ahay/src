function [ D1 ] = fxymssa_win(D,flow,fhigh,dt,N,verb,n1win,n2win,n3win,r1,r2,r3)
%  FXYMSSA_WIN: Local RR method
%
%  IN   D:   	intput 3D data
%       flow:   processing frequency range (lower)
%       fhigh:  processing frequency range (higher)
%       dt:     temporal sampling interval
%       N:      number of singular value to be preserved
%       verb:   verbosity flag (default: 0)
%
%  OUT  D1:  	output data
%
%  Copyright (C) 2013 The University of Texas at Austin
%  Copyright (C) 2013 Yangkang Chen
%  Modified 2015 by Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%  REFERENCES
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
% 
% DEMO: test/test_win3d_fxymssa.m

if nargin==0
    error('Input data must be provided!');
end

[nt,nx,ny]=size(D);

if nargin==1
    flow=1;
    fhigh=124;
    dt=0.004;
    N=1;
    verb=0;
    n1win=nt;
    n2win=nx;
    n3win=ny;
    r1=0.5;
    r2=0.5;
    r3=0.5;
end;

param.dt=dt;
param.flow=flow;
param.fhigh=fhigh;
param.N=N;
param.verb=verb;

D1=win3d(@localfxymssa, param, D, n1win, n2win, n3win, r1, r2, r3);






return



