function Synth(clean,noisy,obs,rr,drr,odrr)
% Author      : Yangkang Chen
% Date        : Oct, 2020
% 
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2019 The University of Texas at Austin
%  Copyright (C) 2019 Yangkang Chen
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
 
%% reconstruct (RR)
flow=1;fhigh=100;dt=0.004;N=10;Niter=10;mode=1;verb=1;iflb=0;a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
NN=Inf;
tic
d1=rr5d_lb_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
toc

%% reconstruct (DRR)
NN=4;
tic
d2=drr5d_lb_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc

%% reconstruct (ODRR)
iflb=2;
tic
d3=drr5d_lb_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc


% figure;
% subplot(4,1,1);imagesc(1:100,[1:100]*4,reshape(d(:,:,:,1,1),nt,10*10));caxis([-0.3,0.3]);colormap(seis);
% xlabel('Trace');ylabel('Time (ms)');title('(a)');
% 
% subplot(4,1,2);imagesc(1:100,[1:100]*4,reshape(d0(:,:,:,1,1),nt,10*10));caxis([-0.3,0.3]);colormap(seis);
% xlabel('Trace');ylabel('Time (ms)');title('(b)');
% 
% subplot(4,1,3);imagesc(1:100,[1:100]*4,reshape(d1(:,:,:,1,1),nt,10*10));caxis([-0.3,0.3]);colormap(seis);
% xlabel('Trace');ylabel('Time (ms)');title('(c)');
% 
% subplot(4,1,4);imagesc(1:100,[1:100]*4,reshape(d2(:,:,:,1,1),nt,10*10));caxis([-0.3,0.3]);colormap(seis);
% xlabel('Trace');ylabel('Time (ms)');title('(d)');
% 
% 
% figure;
% subplot(4,1,1);imagesc(reshape(d(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);colormap(seis);
% xlabel('Trace');ylabel('Time (ms)');title('(a)');
% 
% subplot(4,1,2);imagesc(reshape(d0(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);colormap(seis);
% xlabel('Trace');ylabel('Time (ms)');title('(b)');
% 
% subplot(4,1,3);imagesc(reshape(d1(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);colormap(seis);
% xlabel('Trace');ylabel('Time (ms)');title('(c)');
% 
% subplot(4,1,4);imagesc(reshape(d2(:,5,5,:,:),nt,10*10));caxis([-0.3,0.3]);colormap(seis);
% xlabel('Trace');ylabel('Time (ms)');title('(d)');

yc_snr(d(:),d0(:))
yc_snr(d(:),d1(:))
yc_snr(d(:),d2(:))
yc_snr(d(:),d3(:))

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

rsf_create(odrr,size(d3)');
rsf_write(d3,odrr);






