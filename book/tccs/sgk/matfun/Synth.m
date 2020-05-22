function Synth(clean,noisy,sgk,ksvd)
% Author      : Yangkang Chen
%               Zhejiang University
%         
% Date        : March, 2020

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
% 
% 
% Note that the original synthetic data cannot be found now
% so I use a new one that is often used in my work
% e.g., in 
% Chen et al., 2016, An open-source Matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
% Huang et al., 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
% 
% In addition
% All codes have been significantly rewritten and upgraded so that the general performance of the codes are much improved. Thus the computational time comparison is not exactly the same as previously reported. But, as you will see, SGK is still much faster than the KSVD method without losing much of the effectiveness.
% 
% Key Reference:
% Chen, Y., 2017, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 209, 21-31.
%
% Other related references (e.g., introducing the subroutines)
% Siahsar, M. A. N., Gholtashi, S., Kahoo, A. R., W. Chen, and Y. Chen, 2017, Data-driven multi-task sparse dictionary learning for noise attenuation of 3D seismic data, Geophysics, 82, V385-V396.
% Siahsar, M. A. N., V. Abolghasemi, and Y. Chen, 2017, Simultaneous denoising and interpolation of 2D seismic data using data-driven non-negative dictionary learning, Signal Processing, 141, 309-321.
% Chen, Y., M. Zhang, M. Bai, and W. Chen, 2019, Improving the signal-to-noise ratio of seismological datasets by unsupervised machine learning, Seismological Research Letters, 90, 1552-1564.
% Chen, Y., S. Zu, W. Chen, M. Zhang, and Z. Guan, 2019, Learning the blending spikes using sparse dictionaries, Geophysical Journal International, 218, 1379?1397. 
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51?KS61.
% Chen, Y., S. Fomel, 2015, Random noise attenuation using local signal-and-noise orthogonalization, Geophysics, 80, WD1-WD9.
% Chen, Y., J. Ma, and S. Fomel, 2016, Double-sparsity dictionary for seismic noise attenuation, Geophysics, 81, V17-V30.
% Zu, S., H. Zhou, R. Wu, M. Jiang, and Y. Chen, 2019, Dictionary learning based on dip patch selection training for random noise attenuation, Geophysics, 84, V169?V183.
% Zu, S., H. Zhou, R. Wu, and Y. Chen, 2019, Hybrid-sparsity constrained dictionary learning for iterative deblending of extremely noisy simultaneous-source data, IEEE Transactions on Geoscience and Remote Sensing, 57, 2249-2262.
% etc. 

%% create data
a1=zeros(300,20);
[n,m]=size(a1);
a3=a1;
a4=a1;

k=0;
a=0.1;
b=1;
for t=-0.055:0.002:0.055
    k=k+1;
    b1(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
    b2(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b3(k)=(1-2*(pi*40*t).^2).*exp(-(pi*40*t).^2);
    b4(k)=(1-2*(pi*30*t).^2).*exp(-(pi*30*t).^2);
end
for i=1:m
  t1(i)=round(140);
  t3(i)=round(-6*i+180);
  t4(i)=round(6*i+10);
  a1(t1(i):t1(i)+k-1,i)=b1; 
  a3(t3(i):t3(i)+k-1,i)=b1; 
  a4(t4(i):t4(i)+k-1,i)=b1;
end

temp=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
for j=1:20
    a4=zeros(300,20);
    for i=1:m
  t4(i)=round(6*i+10+3*j); 
  a4(t4(i):t4(i)+k-1,i)=b1;
  
  t1(i)=round(140-2*j);
  a1(t1(i):t1(i)+k-1,i)=b1;
    end
    shot(:,:,j)=a1(1:300,:)+a3(1:300,:)+a4(1:300,:);
end
plane3d=shot;
dc=yc_scale(plane3d,3);
dc=dc(51:225,:,:);
%figure;imagesc(reshape(dc,175,20*20));colormap(seis);

%% adding noise
randn('state',201617);
d=dc+0.1*randn(size(dc));
[n1,n2]=size(d);
%figure;imagesc([dc(:,:,10),d(:,:,10)]);colormap(seis);

%% denoise using SGK
param=struct('T',2,'niter',10,'mode',1,'K',64);
mode=1;l1=4;l2=4;l3=4;s1=2;s2=2;s3=2;perc=1; 
tic
[d1,D,G,D0]=yc_sgk_denoise(d,mode,[l1,l2,l3],[s1,s2,s3],perc,param);
toc
%figure;imagesc([dc(:,:,10),d(:,:,10),d1(:,:,10),d(:,:,10)-d1(:,:,10)]);colormap(seis);

%% benchmark with KSVD
tic
[d2,D2,G2,D02]=yc_ksvd_denoise(d,mode,[l1,l2,l3],[s1,s2,s3],perc,param);
toc
%figure;imagesc([dc(:,:,10),d(:,:,10),d2(:,:,10),d(:,:,10)-d2(:,:,10)]);colormap(seis);

%% SNR
yc_snr(dc,d,2)
yc_snr(dc,d1,2)
yc_snr(dc,d2,2)

rsf_create(clean,size(dc)');
rsf_write(dc,clean);

rsf_create(noisy,size(d)');
rsf_write(d,noisy);

rsf_create(sgk,size(d1)');
rsf_write(d1,sgk);

rsf_create(ksvd,size(d2)');
rsf_write(d2,ksvd);

