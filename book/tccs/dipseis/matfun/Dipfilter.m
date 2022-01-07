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
% parameters definition
 flow=5;
 fhigh=120;
 seed=201314;
 nt=512;dt=0.004;
 nx=256;dx=1;
 f0=40;
 verb=0;
 
%d1  =    linear_events(0.004,40,2,[0:10:10*79],1,0,1,2,2);
d01 =   levents(dt,f0,(nt-1)*0.004,1:nx,0.7,0,0.5,200,2,seed);
% figure;imagesc(d1);
d02 =   levents(dt,f0,(nt-1)*0.004,1:nx,0.8,0.001,0.5,200,2,seed);
% figure;imagesc(d2);
d03 =   levents(dt,f0,(nt-1)*0.004,1:nx,0.3,0.004,0.5,200,2,seed);
d04 =   levents(dt,f0,(nt-1)*0.004,1:nx,0.5,0.005,0.5,200,2,seed);
% figure;imagesc(d3);


d0=d01+d02+d03+d04;d0=d0/max(max(d0));
randn('state',seed);
d1=d0+0.05*randn(nt,nx);
% figure;imagesc([d0,d1]);

dip1=d1-fxemd(d1,5, 120, 0.004, 1, verb);
dip2=d1-fxemd(d1,5, 120, 0.004, 2, verb)-dip1;
dip3=d1-fxemd(d1,5, 120, 0.004, 3, verb)-dip1-dip2;
dip4=d1-fxemd(d1,5, 120, 0.004, 4, verb)-dip1-dip2-dip3;
res=d1-dip1-dip2-dip3-dip4;

dips=[d1,dip1,dip2,dip3,dip4,res];

% from Matlab to Madagascar
rsf_create(plane_clean,size(d0)');
rsf_write(d0,plane_clean);

rsf_create(plane_noisy,size(d1)');
rsf_write(d1,plane_noisy);

rsf_create(planes,size(dips)');
rsf_write(dips,planes);


%figure;imagesc([[d1,dip1,dip2];[dip3,dip4,res]]);
% plane1=dip1+dip2;
% plane2=dip3+dip4+res;
% figure;imagesc([plane1,plane2]);

%figure;imagesc([plane1,plane2,plane3]);

%figure;imagesc([d1,dip1,dip2;[dip3,dip4,res]]);
