function Synth(clean,noisy,obs,rr,drr)
% Author      : Yangkang Chen
% Date        : Feb, 2020
% 
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2016 The University of Texas at Austin
%  Copyright (C) 2016 Yangkang Chen
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

% RR
flow=0;fhigh=100;dt=0.004;N=6;Niter=10;mode=1;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=rr5d_lb_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);


% DRR
NN=2; %or N=3;NN=3; the SNR is even higher
d2=drr5d_lb_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);

yc_snr(d(:),d0(:))
yc_snr(d(:),d1(:))
yc_snr(d(:),d2(:))

%% save data
% fid1=fopen('synth_clean5d.bin','w');
% fwrite(fid1,d,'float');
% fid2=fopen('synth_noisy5d.bin','w');
% fwrite(fid2,dn,'float');
% fid3=fopen('synth_obs5d.bin','w');
% fwrite(fid3,d0,'float');
% fid4=fopen('synth_rr5d.bin','w');
% fwrite(fid4,d1,'float');
% fid5=fopen('synth_irr5d.bin','w');
% fwrite(fid5,d2,'float');


%% from Matlab to Madagascar
rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(obs,size(d0)');
rsf_write(d0,obs);

rsf_create(rr,size(d1)');
rsf_write(d1,rr);

rsf_create(drr,size(d2)');
rsf_write(d2,drr);






