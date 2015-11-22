 function field(data,fxmssa,fxemd,fxemdmssa)
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

% parameters definition
 

% allocate memory
 D = zeros(1024,471);
 
 % from Madagascar to Matlab
 rsf_read(D,data);

 flow=5;
 fhigh=245;
 dt=0.002;
 N1=3;
 N2=10;
 verb=1;

 %% Using fx_mssa
  d11= fx_mssa(D,flow,fhigh,dt,N2,verb);  			
 
 % processing in a separated window 
 % i# vertical i## horizontal ABCD (anti-clockwise) 
 % w1 vertical w2 horizontal
 i1=340; i11=300;  i3=340; i33=420;
 i2=440; i22=300;  i4=440; i44=420;
 w1=20;
 w2=25;
 
 d1t=fx_mssa([D(i1:i2,i11:i33)-d11(i1:i2,i11:i33)],flow,fhigh,dt,5,verb)+d11(i1:i2,i11:i33); % MSSA on difference section.
 d1=d11;
 d1(i1+w2:i2-w2+1,i11+w1:i33-w1+1)=d1t(w2+1:i2-i1+2-w2,w1+1:i33-i11+2-w1);
 d1(i1:i1+w2-1,i11+w1:i33-w1)=d11(i1:i1+w2-1,i11+w1:i33-w1).*([w2:-1:1]'./w2*ones(1,i33-i22-2*w1+1))+d1t(1:w2,w1+1:i33-i22-w1+1).*([1:1:w2]'./w2*ones(1,i33-i22-2*w1+1));
 d1(i2-w2+1:i2,i11+w1:i33-w1)=d11(i2-w2+1:i2,i11+w1:i33-w1).*([1:1:w2]'./w2*ones(1,i33-i22-2*w1+1))+d1t(i2-i1+2-w2:i2-i1+1,w1+1:i33-i22-w1+1).*([w2:-1:1]'/w2*ones(1,i33-i22-2*w1+1));
 d1(i1:i2,i11:i11+w1-1)=d11(i1:i2,i11:i11+w1-1).*(ones(i2-i1+1,1)*[w1:-1:1]/w1)+d1t(:,1:w1).*(ones(i2-i1+1,1)*[1:1:w1]/w1);
 d1(i1:i2,i33-w1+1:i33)=d11(i1:i2,i33-w1+1:i33).*(ones(i2-i1+1,1)*[1:1:w1]/w1)+d1t(:,i33-i22+2-w1:i33-i22+1).*(ones(i2-i1+1,1)*[w1:-1:1]/w1);
 %figure; imagesc([D,d1,D-d1],[-0.5,0.5]);colormap(seismic);

 %% Using fx_emd 
 d2= fx_emd(D,flow,fhigh,dt,N1,verb);  			%using fx_emd
 %figure; imagesc([D,d2,D-d2],[-0.5,0.5]);colormap(seismic);
 

 %% Using fx_emdmssa
 d31=d2;
 d3t=fx_mssa([D(i1:i2,i11:i33)-d31(i1:i2,i11:i33)],flow,fhigh,dt,5,verb)+d31(i1:i2,i11:i33); % MSSA on difference section.
 d3=d31;
 d3(i1+w2:i2-w2+1,i11+w1:i33-w1+1)=d3t(w2+1:i2-i1+2-w2,w1+1:i33-i11+2-w1);
 d3(i1:i1+w2-1,i11+w1:i33-w1)=d31(i1:i1+w2-1,i11+w1:i33-w1).*([w2:-1:1]'./w2*ones(1,i33-i22-2*w1+1))+d3t(1:w2,w1+1:i33-i22-w1+1).*([1:1:w2]'./w2*ones(1,i33-i22-2*w1+1));
 d3(i2-w2+1:i2,i11+w1:i33-w1)=d31(i2-w2+1:i2,i11+w1:i33-w1).*([1:1:w2]'./w2*ones(1,i33-i22-2*w1+1))+d3t(i2-i1+2-w2:i2-i1+1,w1+1:i33-i22-w1+1).*([w2:-1:1]'/w2*ones(1,i33-i22-2*w1+1));
 d3(i1:i2,i11:i11+w1-1)=d31(i1:i2,i11:i11+w1-1).*(ones(i2-i1+1,1)*[w1:-1:1]/w1)+d3t(:,1:w1).*(ones(i2-i1+1,1)*[1:1:w1]/w1);
 d3(i1:i2,i33-w1+1:i33)=d31(i1:i2,i33-w1+1:i33).*(ones(i2-i1+1,1)*[1:1:w1]/w1)+d3t(:,i33-i22+2-w1:i33-i22+1).*(ones(i2-i1+1,1)*[w1:-1:1]/w1);
 
 %figure; imagesc([D,d3,D-d3],[-0.5,0.5]);colormap(seismic);
 
 rsf_create(fxmssa,size(d1)');
 rsf_write(d1,fxmssa);
 
 rsf_create(fxemd,size(d2)');
 rsf_write(d2,fxemd);
 
 rsf_create(fxemdmssa,size(d3)');
 rsf_write(d3,fxemdmssa);
 
 
 
 
 
 
 
