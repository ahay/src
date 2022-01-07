function Hankel(H_clean,H_obs,H_rr,H_drr,H_odrr,H_rr10,H_drr10,H_odrr10)
% Author      : Yangkang Chen
% Date        : Oct, 2020
% 
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2020 The University of Texas at Austin
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
% 
%  Reference:   Chen et al., GJI, 2016; Chen et al., GJI, 2020


%% load data
load yc_synth5d.mat
d=data5d;d=d/max(max(max(max(max(d)))));
[nt,nhx,nhy,nx,ny]=size(d);

%% simultaneous denoising and reconstruction
randn('state',201314);
var=0.25;
dn=d+var*randn(size(d));

%% decimate
[nt,nhx,nhy,nx,ny]=size(d);
ratio=0.3;
mask=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
d0=dn.*mask;

%% reconstruct
flow=0;fhigh=100;dt=0.004;N=10;NN=4;Niter=10;mode=1;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
tic
[d1,M_c,M_obs,M_rr,M_rr10]=drr5d_lb_recon_H(d0,d,mask,flow,fhigh,30,dt,N,Inf,Niter,eps,verb,mode,iflb,a);
toc

tic
[d2,M_c,M_obs,M_drr,M_drr10]=drr5d_lb_recon_H(d0,d,mask,flow,fhigh,30,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc

tic
[d3,M_c,M_obs,M_odrr,M_odrr10]=drr5d_lb_recon_H(d0,d,mask,flow,fhigh,30,dt,N,NN,Niter,eps,verb,mode,2,a);
toc

%% from Matlab to Madagascar
rsf_create(H_clean,size(M_c)');
rsf_write(abs(M_c),H_clean);

rsf_create(H_obs,size(M_obs)');
rsf_write(abs(M_obs),H_obs);

rsf_create(H_rr,size(M_rr)');
rsf_write(abs(M_rr),H_rr);

rsf_create(H_drr,size(M_drr)');
rsf_write(abs(M_drr),H_drr);

rsf_create(H_odrr,size(M_odrr)');
rsf_write(abs(M_odrr),H_odrr);

rsf_create(H_rr10,size(M_rr10)');
rsf_write(abs(M_rr10),H_rr10);

rsf_create(H_drr10,size(M_drr10)');
rsf_write(abs(M_drr10),H_drr10);

rsf_create(H_odrr10,size(M_odrr10)');
rsf_write(abs(M_odrr10),H_odrr10);

