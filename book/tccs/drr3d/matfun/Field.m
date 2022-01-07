function Field(obs,mssa,dmssa)
% Author      : Yangkang Chen
%               Zhejiang University
%         
% Date        : Feb, 2020

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Yangkang Chen
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
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  

fid=fopen('real3d.bin','r');
d=fread(fid,[300,1000],'float');
d=reshape(d,300,100,10);

%% without noise
dn=d;

%% decimate
[nt,nx,ny]=size(d);
ratio=0.5;
mask=genmask(reshape(d,nt,nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny);
d0=dn.*mask;

%% simultaneous denoising and reconstruction
% adding noise
randn('state',201314);
var=0.1;
dn=d+var*randn(size(d));
d0=dn.*mask;

flow=0;fhigh=125;dt=0.004;N=10;Niter=10;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=fxymssa_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);

% figure;imagesc([d0(:,:,5),d1(:,:,5)]);
%% from Matlab to Madagascar
rsf_create(obs,size(d0)');
rsf_write(d0,obs);

rsf_create(mssa,size(d1)');
rsf_write(d1,mssa);

flow=0;fhigh=125;dt=0.004;N=10;Niter=10;mode=1;verb=1;NN=3;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d2=fxydmssa_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,a);

rsf_create(dmssa,size(d2)');
rsf_write(d2,dmssa);



