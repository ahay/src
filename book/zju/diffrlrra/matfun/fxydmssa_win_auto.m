function [ D1 ] = fxydmssa_win_auto(D,flow,fhigh,dt,N,NN,verb,n1win,n2win,n3win,r1,r2,r3,mode)
%  FXY_DMSSA: F-XY domain localized damped multichannel singular spectrum analysis (MSSA)
%  with automatic rank selection
%
%  IN   D:   	intput 3D data
%       flow:   processing frequency range (lower)
%       fhigh:  processing frequency range (higher)
%       dt:     temporal sampling interval
%       N:      number of singular value to be preserved
%       NN:     damping factor
%       verb:   verbosity flag (default: 0)
%
%  OUT  D1:  	output data
%
%  Copyright (C) 2013 The University of Texas at Austin
%  Copyright (C) 2013 Yangkang Chen
%  Modified 2015 by Yangkang Chen
%  Modified 2020 by Yangkang Chen
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
%  Reference:   Simultaneous seismic data denoising and reconstruction via multichannel
%               singular spectrum analysis, Geophysics, 2011, 76, V25-V32
%
% DEMO: test/test_win3d_fxydmssa_autorank.m %or so

if nargin==0
    error('Input data must be provided!');
end

[nt,nx,ny]=size(D);

if nargin==1
    flow=1;
    fhigh=124;
    dt=0.004;
    N=1;
    N=2;
    verb=0;
    n1win=nt;
    n2win=nx;
    n3win=ny;
    r1=0.5;
    r2=0.5;
    r3=0.5;
    mode=2;
end;

param.dt=dt;
param.flow=flow;
param.fhigh=fhigh;
param.N=N;
param.NN=NN;
param.verb=verb;
param.mode=mode;

D1=win3d(@localfxydmssa_auto, param, D, n1win, n2win, n3win, r1, r2, r3);






return



