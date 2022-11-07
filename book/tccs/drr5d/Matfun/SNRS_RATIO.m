function SNRS_RATIO(synth_drr_rs)
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
var=0.25;
randn('state',201314);
dn=d+var*randn(size(d));

%% parameters
flow=1;fhigh=80;dt=0.004;N=6;NN=2;Niter=10;mode=1;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing

%% case 1: ratio=0.1 (90% missing)
r1=0.1;
mask1=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),r1,'c',201415);
mask1=reshape(mask1,nt,nhx,nhy,nx,ny);
% %% decimate
d1_0=dn.*mask1;
% %%main processing
tic
d1_1=rr5d_lb_recon(d1_0,mask1,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d1_2=drr5d_lb_recon(d1_0,mask1,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr1_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d1_1,nt,nhx*nhy*nx*ny))
snr1_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d1_2,nt,nhx*nhy*nx*ny))
fprintf('Case 1 finished !!!\n\n');

%% case 2: ratio=0.2 (80% missing)
r2=0.2;
mask2=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),r2,'c',201415);
mask2=reshape(mask2,nt,nhx,nhy,nx,ny);
% %% decimate
d2_0=dn.*mask2;
% %%main processing
tic
d2_1=rr5d_lb_recon(d2_0,mask2,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d2_2=drr5d_lb_recon(d2_0,mask2,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr2_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d2_1,nt,nhx*nhy*nx*ny))
snr2_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d2_2,nt,nhx*nhy*nx*ny))
fprintf('Case 2 finished !!!\n\n');

%% case 3: ratio=0.3 (70% missing)
r3=0.3;
mask3=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),r3,'c',201415);
mask3=reshape(mask3,nt,nhx,nhy,nx,ny);
% %% decimate
d3_0=dn.*mask3;
% %%main processing
tic
d3_1=rr5d_lb_recon(d3_0,mask3,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d3_2=drr5d_lb_recon(d3_0,mask3,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr3_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d3_1,nt,nhx*nhy*nx*ny))
snr3_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d3_2,nt,nhx*nhy*nx*ny))
fprintf('Case 3 finished !!!\n\n');

%% case 4: ratio=0.4 (60% missing)
r4=0.4;
mask4=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),r4,'c',201415);
mask4=reshape(mask4,nt,nhx,nhy,nx,ny);
% %% decimate
d4_0=dn.*mask4;
% %%main processing
tic
d4_1=rr5d_lb_recon(d4_0,mask4,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d4_2=drr5d_lb_recon(d4_0,mask4,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr4_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d4_1,nt,nhx*nhy*nx*ny))
snr4_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d4_2,nt,nhx*nhy*nx*ny))
fprintf('Case 4 finished !!!\n\n');

%% case 5: ratio=0.5 (50% missing)
r5=0.5;
mask5=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),r5,'c',201415);
mask5=reshape(mask5,nt,nhx,nhy,nx,ny);
% %% decimate
d5_0=dn.*mask5;
% %%main processing
tic
d5_1=rr5d_lb_recon(d5_0,mask5,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d5_2=drr5d_lb_recon(d5_0,mask5,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr5_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d5_1,nt,nhx*nhy*nx*ny));
snr5_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d5_2,nt,nhx*nhy*nx*ny));
fprintf('Case 5 finished !!!\n\n');

%% case 6: ratio=0.6 (40% missing)
r6=0.6;
mask6=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),r6,'c',201415);
mask6=reshape(mask6,nt,nhx,nhy,nx,ny);
% %% decimate
d6_0=dn.*mask6;
% %%main processing
tic
d6_1=rr5d_lb_recon(d6_0,mask6,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d6_2=drr5d_lb_recon(d6_0,mask6,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr6_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d6_1,nt,nhx*nhy*nx*ny));
snr6_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d6_2,nt,nhx*nhy*nx*ny));
fprintf('Case 6 finished !!!\n\n');

%% case 7: ratio=0.7 (30% missing)
r7=0.7;
mask7=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),r7,'c',201415);
mask7=reshape(mask7,nt,nhx,nhy,nx,ny);
% %% decimate
d7_0=dn.*mask7;
% %%main processing
tic
d7_1=rr5d_lb_recon(d7_0,mask7,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d7_2=drr5d_lb_recon(d7_0,mask7,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr7_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d7_1,nt,nhx*nhy*nx*ny));
snr7_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d7_2,nt,nhx*nhy*nx*ny));
fprintf('Case 7 finished !!!\n\n');

%% case 6: ratio=0.6 (40% missing)
r8=0.8;
mask8=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),r8,'c',201415);
mask8=reshape(mask8,nt,nhx,nhy,nx,ny);
% %% decimate
d8_0=dn.*mask8;
% %%main processing
tic
d8_1=rr5d_lb_recon(d8_0,mask8,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d8_2=drr5d_lb_recon(d8_0,mask8,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr8_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d8_1,nt,nhx*nhy*nx*ny));
snr8_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d8_2,nt,nhx*nhy*nx*ny));
fprintf('Case 8 finished !!!\n\n');

%% plot SNR
rs=[r1 r2 r3 r4 r5 r6 r7 r8];
snrs1=[snr1_1 snr2_1 snr3_1 snr4_1 snr5_1 snr6_1 snr7_1 snr8_1];
snrs2=[snr1_2 snr2_2 snr3_2 snr4_2 snr5_2 snr6_2 snr7_2 snr8_2];
%figure;plot(rs,snrs1,'ro');hold on;
%plot(rs,snrs2,'bv');

%% from Matlab to Madagascar
snrs=[rs',snrs1',snrs2'];
rsf_create(synth_drr_rs,size(snrs)');
rsf_write(snrs,synth_drr_rs);

