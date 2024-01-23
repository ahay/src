function Hyper(data,clean,noisy,fk,rr,drr,odrr)
% Author      : Min Bai and Yangkang Chen
%               Zhejiang University
% 
% Date        : Feb, 2020
% 
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Min Bai and Yangkang Chen
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

%% generate synthetic data


shot3d=zeros(126,128*128);

rsf_read(shot3d,data);
shot3d=reshape(shot3d,126,128,128);
shot3d=shot3d(:,65:96,65:96);
d=yc_scale(shot3d,3);


%% adding noise
randn('state',201314);
var=0.2;
dn=d+var*randn(size(d));
yc_snr(d,dn,2)
%snr:-8.39

%% denoise
flow=0;fhigh=80;dt=0.004;N=20;verb=0;

%% RR
mode=1;
tic
d1=fxyodrr(dn(:,:,:),flow,fhigh,dt,N,mode,verb,6);
toc
%figure;imagesc([d(:,:,9),dn(:,:,9),d1(:,:,9),dn(:,:,9)-d1(:,:,9)]);
%colormap(seis);
yc_snr(d,d1,2)
%about 1.69s
%snr:8.27

%% DRR
mode=2;
tic
d2=fxyodrr(dn(:,:,:),flow,fhigh,dt,N,mode,verb,6);
toc
%figure;imagesc([d(:,:,9),dn(:,:,9),d4(:,:,9),dn(:,:,9)-d4(:,:,9)]);
%colormap(seis);
yc_snr(d,d2,2)
%about 1.71s
%SNR: 9.58

%% ODRR
mode=3;
tic
d3=fxyodrr(dn(:,:,:),flow,fhigh,dt,N,mode,verb,6);
toc
yc_snr(d,d3,2)
%about 3.34s
%SNR: 9.65

%% FK
tic
d4=yc_fkt(dn,'ps',40);
toc
yc_snr(d,d4,2)
%about 0.015s
%SNR: 7.05


%% from Matlab to Madagascar
rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(rr,size(d1)');
rsf_write(d1,rr);

rsf_create(drr,size(d2)');
rsf_write(d2,drr);

rsf_create(odrr,size(d3)');
rsf_write(d3,odrr);

rsf_create(fk,size(d4)');
rsf_write(d4,fk);



