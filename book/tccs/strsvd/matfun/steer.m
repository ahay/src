function steer(dn,dh,dhc,dc,s)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Jul, 2014

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
%  
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% SVD after dip-steering
% Reference: 
% Bekara, M., and M. van der Baan, 2007, Local singular value decomposition for signal
%   enhancement of seismic data: Geophysics, 72, V59–V65.


%clear;clc;close all;
  dt = 2./1000;
  tmax = 1.;
  h = [0:10:10*(50-1)];
  tau = [0.1,0.2,0.3,0.6];
  p = [0.0004,-0.0001,0.0001,-0.0001];
  amp = [1.2,-1.,1.,1];
  f0 = 40;
  snr = 5;
  L = 1;
  seed=201314;
  
  %d=linear_events(dt,f0,tmax,h,tau,p,amp,snr,L,seed);
  d=linear_events(dt,f0,tmax,h,0.4,0.0006,1,snr,L,seed);
  %figure;imagesc(d);
  
  %% dip steering
  ref=d(:,1);
  shift=dip_steering(d,ref);
  d1=seisdither(d,shift);
  
  %% inverse dip steering
  d2=seisdither(d1,-shift);
  
  %% SVD after dip steering
  [u,e,v]=svd(d1);
  r=1;
  d1_svd=u(:,1:r)*e(1:r,1:r)*v(:,1:r)';
  d2_svd=seisdither(d1_svd,-shift);

% from Matlab to Madagascar
dnoisy='dd.rsf';
rsf_create(dn,size(d)');
rsf_write(d,dn);

rsf_create(dh,size(d1)');
rsf_write(d1,dh);

rsf_create(dhc,size(d1_svd)');
rsf_write(d1_svd,dhc);

rsf_create(dc,size(d2_svd)');
rsf_write(d2_svd,dc);
 
rsf_create(s,size(shift)');
rsf_write(shift,s); 
