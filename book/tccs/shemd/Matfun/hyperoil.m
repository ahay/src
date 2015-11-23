 function synth(clean,noisy,fxmssa,fxemd,fxemdmssa)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Nov, 2013

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2013 The University of Texas at Austin
%  Copyright (C) 2013 Yangkang Chen
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

% create synthetic data
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
  
  d=hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L,seed);
d=d/max(max(d));
  
D=d;
%rng(201313);
randn('state',201313);
D=D+0.1*randn(size(d));
%figure;wigb(D);

% parameters definition
 flow=5;
 fhigh=245;
 dt=0.002;
 N1=1;
 N2=6;
 verb=1;

 %% Using fx_mssa
 d1= fx_mssa(D,flow,fhigh,dt,N2,verb);  			
 %figure; imagesc([D,d1,D-d1],[-0.5,0.5]);colormap(seismic);

 %% Using fx_emd 
 d2= fx_emd(D,flow,fhigh,dt,N1,verb);  			
 %figure; imagesc([D,d2,D-d2],[-0.5,0.5]);colormap(seismic);
 
 %% Using fx_emdmssa
 d3=fx_emdmssa(D,flow,fhigh,dt,N1,1,verb);
 %figure; imagesc([D,d3,D-d3],[-0.5,0.5]);colormap(seismic);
 
 %% Using fx_emdmssa
 d31=d2;
 % processing in a separated window 1
 i1=50; i11=1;  i3=50; i33=25;
 i2=170; i22=1;  i4=170; i44=25;
 w1=2;
 w2=2;

 d3t=fx_mssa([D(i1:i2,i11:i33)-d31(i1:i2,i11:i33)],flow,fhigh,dt,2,verb)+d31(i1:i2,i11:i33); % MSSA on difference section.
 d3=d31;
 d3(i1+w2:i2-w2+1,i11+w1:i33-w1+1)=d3t(w2+1:i2-i1+2-w2,w1+1:i33-i11+2-w1);
 d3(i1:i1+w2-1,i11+w1:i33-w1)=d31(i1:i1+w2-1,i11+w1:i33-w1).*([w2:-1:1]'./w2*ones(1,i33-i22-2*w1+1))+d3t(1:w2,w1+1:i33-i22-w1+1).*([1:1:w2]'./w2*ones(1,i33-i22-2*w1+1));
 d3(i2-w2+1:i2,i11+w1:i33-w1)=d31(i2-w2+1:i2,i11+w1:i33-w1).*([1:1:w2]'./w2*ones(1,i33-i22-2*w1+1))+d3t(i2-i1+2-w2:i2-i1+1,w1+1:i33-i22-w1+1).*([w2:-1:1]'/w2*ones(1,i33-i22-2*w1+1));
 d3(i1:i2,i11:i11+w1-1)=d31(i1:i2,i11:i11+w1-1).*(ones(i2-i1+1,1)*[w1:-1:1]/w1)+d3t(:,1:w1).*(ones(i2-i1+1,1)*[1:1:w1]/w1);
 d3(i1:i2,i33-w1+1:i33)=d31(i1:i2,i33-w1+1:i33).*(ones(i2-i1+1,1)*[1:1:w1]/w1)+d3t(:,i33-i22+2-w1:i33-i22+1).*(ones(i2-i1+1,1)*[w1:-1:1]/w1);


 % processing in a separated window 2
 i1=50; i11=25;  i3=50; i33=76;
 %i2=155; i22=78;  i4=155; i44=99;
 i2=350; i22=25;  i4=350; i44=76;
 w1=5;
 w2=5;

 d3tt=fx_mssa([D(i1:i2,i11:i33)-d31(i1:i2,i11:i33)],flow,fhigh,dt,3,verb)+d31(i1:i2,i11:i33); % MSSA on difference section.
 d3(i1+w2:i2-w2+1,i11+w1:i33-w1+1)=d3tt(w2+1:i2-i1+2-w2,w1+1:i33-i11+2-w1);
 d3(i1:i1+w2-1,i11+w1:i33-w1)=d31(i1:i1+w2-1,i11+w1:i33-w1).*([w2:-1:1]'./w2*ones(1,i33-i22-2*w1+1))+d3tt(1:w2,w1+1:i33-i22-w1+1).*([1:1:w2]'./w2*ones(1,i33-i22-2*w1+1));
 d3(i2-w2+1:i2,i11+w1:i33-w1)=d31(i2-w2+1:i2,i11+w1:i33-w1).*([1:1:w2]'./w2*ones(1,i33-i22-2*w1+1))+d3tt(i2-i1+2-w2:i2-i1+1,w1+1:i33-i22-w1+1).*([w2:-1:1]'/w2*ones(1,i33-i22-2*w1+1));
 d3(i1:i2,i11:i11+w1-1)=d31(i1:i2,i11:i11+w1-1).*(ones(i2-i1+1,1)*[w1:-1:1]/w1)+d3tt(:,1:w1).*(ones(i2-i1+1,1)*[1:1:w1]/w1);
 d3(i1:i2,i33-w1+1:i33)=d31(i1:i2,i33-w1+1:i33).*(ones(i2-i1+1,1)*[1:1:w1]/w1)+d3tt(:,i33-i22+2-w1:i33-i22+1).*(ones(i2-i1+1,1)*[w1:-1:1]/w1);

 
 rsf_create(clean,size(d)');
 rsf_write(d,clean);

 rsf_create(noisy,size(D)');
 rsf_write(D,noisy);

 rsf_create(fxmssa,size(d1)');
 rsf_write(d1,fxmssa);
 
 rsf_create(fxemd,size(d2)');
 rsf_write(d2,fxemd);
 
 rsf_create(fxemdmssa,size(d3)');
 rsf_write(d3,fxemdmssa);
 
 
 
 
 
 
 
