function SNRS_SIGMA(synth_drr_snrs)
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
ratio=0.3;
mask=yc_genmask(reshape(d,nt,nhx*nhy*nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nhx,nhy,nx,ny);
randn('state',201314);

%% parameters
flow=1;fhigh=80;dt=0.004;N=6;NN=2;Niter=10;mode=1;verb=1;iflb=0;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing

%% case 1: sigma=0.05
var1=0.05;
d1_n=d+var1*randn(size(d));
% decimate
d1_0=d1_n.*mask;
% main processing
tic
d1_1=rr5d_lb_recon(d1_0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d1_2=drr5d_lb_recon(d1_0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr1_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d1_1,nt,nhx*nhy*nx*ny));
snr1_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d1_2,nt,nhx*nhy*nx*ny));
fprintf('Case 1 finished !!!\n\n');

%% case 2: sigma=0.1
var2=0.1;
d2_n=d+var2*randn(size(d));
% decimate
d2_0=d2_n.*mask;
% main processing
tic
d2_1=rr5d_lb_recon(d2_0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d2_2=drr5d_lb_recon(d2_0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr2_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d2_1,nt,nhx*nhy*nx*ny));
snr2_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d2_2,nt,nhx*nhy*nx*ny));
fprintf('Case 2 finished !!!\n\n');

%% case 3: sigma=0.25
var3=0.25;
d3_n=d+var3*randn(size(d));
% decimate
d3_0=d3_n.*mask;
% main processing
tic
d3_1=rr5d_lb_recon(d3_0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d3_2=drr5d_lb_recon(d3_0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr3_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d3_1,nt,nhx*nhy*nx*ny));
snr3_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d3_2,nt,nhx*nhy*nx*ny));
fprintf('Case 3 finished !!!\n\n');

%% case 4: sigma=0.5
var4=0.5;
d4_n=d+var4*randn(size(d));
% decimate
d4_0=d4_n.*mask;
% main processing
tic
d4_1=rr5d_lb_recon(d4_0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d4_2=drr5d_lb_recon(d4_0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr4_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d4_1,nt,nhx*nhy*nx*ny));
snr4_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d4_2,nt,nhx*nhy*nx*ny));
fprintf('Case 4 finished !!!\n\n');

%% case 5: sigma=0.75
var5=0.75;
d5_n=d+var5*randn(size(d));
% decimate
d5_0=d5_n.*mask;
% main processing
tic
d5_1=rr5d_lb_recon(d5_0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d5_2=drr5d_lb_recon(d5_0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr5_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d5_1,nt,nhx*nhy*nx*ny));
snr5_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d5_2,nt,nhx*nhy*nx*ny));
fprintf('Case 5 finished !!!\n\n');

%% case 6: sigma=1
var6=1;
d6_n=d+var6*randn(size(d));
% decimate
d6_0=d6_n.*mask;
% main processing
tic
d6_1=rr5d_lb_recon(d6_0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,iflb,a);
d6_2=drr5d_lb_recon(d6_0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,iflb,a);
toc
snr6_1=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d6_1,nt,nhx*nhy*nx*ny));
snr6_2=yc_snr(reshape(d,nt,nhx*nhy*nx*ny),reshape(d6_2,nt,nhx*nhy*nx*ny));
fprintf('Case 6 finished !!!\n\n');

%% plot SNR
vars=[var1 var2 var3 var4 var5 var6];
snrs1=[snr1_1 snr2_1 snr3_1 snr4_1 snr5_1 snr6_1];
snrs2=[snr1_2 snr2_2 snr3_2 snr4_2 snr5_2 snr6_2];
%figure;plot(vars,snrs1,'ro');hold on;
%plot(vars,snrs2,'bv');

%% from Matlab to Madagascar
snrs=[vars',snrs1',snrs2'];
rsf_create(synth_drr_snrs,size(snrs)');
rsf_write(snrs,synth_drr_snrs);


