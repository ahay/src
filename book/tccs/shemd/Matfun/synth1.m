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
[d1] = linear_events(0.002,40,1.0, 1:50 ,[0.3,0.6],...
    [0,0],[-1.,1],200,5,2013); 


[d2] = linear_events(0.002,40,1.0, 1:50, [0.1],...
    [0.014],[2.],200,5,2013);
d2(1:120,:)=zeros(120,50);
d2(350:501,:)=zeros(152,50);

[d3] = linear_events(0.002,40,1.0,1:50,0.2,0.000,1,200,5,2013);

d3(:,1:15)=zeros(501,15);

[d4] = linear_events(0.002,40,1.0,1:50,0.8,0.000,1,200,5,2013);
d4(:,41:50)=zeros(501,10);

d=d1+d2+d3+d4;
d=d/max(max(d));
D=d;
%rng(201313);
randn('state',201313);
D=D+0.15*randn(size(d));
%figure;wigb(D);

% parameters definition
 flow=5;
 fhigh=245;
 dt=0.002;
 N1=3;
 N2=4;
 verb=1;

 %% Using fx_mssa
 d1= fx_mssa(D,flow,fhigh,dt,N2,verb);  			
 %figure; imagesc([D,d1,D-d1],[-0.5,0.5]);colormap(seismic);

 %% Using fx_emd 
 d2= fx_emd(D,flow,fhigh,dt,N1,verb);  			
 %figure; imagesc([D,d2,D-d2],[-0.5,0.5]);colormap(seismic);
 
 %% Using fx_emdmssa
 %d3=fx_emdmssa(D,flow,fhigh,dt,N1,1,verb);
 %figure; imagesc([D,d3,D-d3],[-0.5,0.5]);colormap(seismic);
 
 %% Using fx_emdmssa
 d31=d2;
 % processing in a separated window 1
 i1=120; i11=8;  i3=120; i33=44;
 i2=350; i22=8;  i4=350; i44=44;
 w1=5;
 w2=2;

 d3t=fx_mssa([D(i1:i2,i11:i33)-d31(i1:i2,i11:i33)],flow,fhigh,dt,1,verb)+d31(i1:i2,i11:i33); % MSSA on difference section.
 d3=d31;
 d3(i1+w2:i2-w2+1,i11+w1:i33-w1+1)=d3t(w2+1:i2-i1+2-w2,w1+1:i33-i11+2-w1);
 d3(i1:i1+w2-1,i11+w1:i33-w1)=d31(i1:i1+w2-1,i11+w1:i33-w1).*([w2:-1:1]'./w2*ones(1,i33-i22-2*w1+1))+d3t(1:w2,w1+1:i33-i22-w1+1).*([1:1:w2]'./w2*ones(1,i33-i22-2*w1+1));
 d3(i2-w2+1:i2,i11+w1:i33-w1)=d31(i2-w2+1:i2,i11+w1:i33-w1).*([1:1:w2]'./w2*ones(1,i33-i22-2*w1+1))+d3t(i2-i1+2-w2:i2-i1+1,w1+1:i33-i22-w1+1).*([w2:-1:1]'/w2*ones(1,i33-i22-2*w1+1));
 d3(i1:i2,i11:i11+w1-1)=d31(i1:i2,i11:i11+w1-1).*(ones(i2-i1+1,1)*[w1:-1:1]/w1)+d3t(:,1:w1).*(ones(i2-i1+1,1)*[1:1:w1]/w1);
 d3(i1:i2,i33-w1+1:i33)=d31(i1:i2,i33-w1+1:i33).*(ones(i2-i1+1,1)*[1:1:w1]/w1)+d3t(:,i33-i22+2-w1:i33-i22+1).*(ones(i2-i1+1,1)*[w1:-1:1]/w1);


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
 
 
 
 
 
 
 
