function Dipfilter(plane_clean,plane_noisy,planes)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Apr, 2014

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


% clear;clc;close all;


%% make model 
 dt = 2./1000;
  tmax = 1.0;
  h = [-500:20:1000];
  tau = [0.1,.4,0.7];
  v = [1500,2400,2300];
  amp = [1., -1.,1];
  f0 = 20;
  snr = 200;
  L = 9;
  seed=2013;
  
  d=hevents(dt,f0,tmax,h,tau,v,amp,snr,L,seed);
d0=d/max(max(d));
[nt,nx]=size(d0);
randn('state',seed);
d1n=d0+0.1*randn(nt,nx);
% figure;imagesc([d0,d1]);

d1=yc_fdct_wrapping_thr(d1n,8);
verb=1;
dip1=d1-fxemd(d1,5, 120, 0.004, 1, verb);
res=d1-dip1;

dips=[d1,dip1,res];

% from Matlab to Madagascar
rsf_create(plane_clean,size(d0)');
rsf_write(d0,plane_clean);

rsf_create(plane_noisy,size(d1n)');
rsf_write(d1n,plane_noisy);

rsf_create(planes,size(dips)');
rsf_write(dips,planes);


%figure;imagesc([[d1,dip1,dip2];[dip3,dip4,res]]);
% plane1=dip1+dip2;
% plane2=dip3+dip4+res;
% figure;imagesc([plane1,plane2]);

%figure;imagesc([plane1,plane2,plane3]);

%figure;imagesc([d1,dip1,dip2;[dip3,dip4,res]]);
