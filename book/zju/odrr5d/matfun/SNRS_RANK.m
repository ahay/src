function SNRS_RANK(synth_drr_ranks)
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
n=randn(size(d));


%% decimate
ranks=4:1:15;

snrs0=zeros(length(ranks),1);
snrs1=zeros(length(ranks),1);
snrs2=zeros(length(ranks),1);
snrs3=zeros(length(ranks),1);

ratio=0.3;
var=0.25;
dn=d+var*n;
mask=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
d0=dn.*mask;

fprintf('N_ranks is %d\n',length(ranks));
for i=1:length(ranks)

N=ranks(i);

%% parameters
flow=5;fhigh=100;dt=0.004;Niter=10;mode=1;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing

%% reconstruct (RR)
NN=Inf;
d1=drr5d_lb_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);

%% reconstruct (DRR)
NN=4;
d2=drr5d_lb_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);

%% reconstruct (ODRR)
iflb=2;
d3=drr5d_lb_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);

snrs0(i)=yc_snr(d(:),d0(:));
snrs1(i)=yc_snr(d(:),d1(:));
snrs2(i)=yc_snr(d(:),d2(:));
snrs3(i)=yc_snr(d(:),d3(:));

fprintf('rank=%g is done\n',ranks(i));
end

%% from Matlab to Madagascar
snrs=[ranks(:),snrs1,snrs2,snrs3];
rsf_create(synth_drr_ranks,size(snrs)');
rsf_write(snrs,synth_drr_ranks);


