function Syn(clean,noisy,fk,rr,drr,odrr)
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
d=plane3d/max(max(max(plane3d)));

%% adding noise
randn('state',20131415);
var=0.2;
dn=d+var*randn(size(d));
fprintf('SNR of Noisy is %g\n',yc_snr(d,dn,2));
%snr:-8.39

%% denoise
flow=0;fhigh=250;dt=0.002;N=3;verb=0;

%% RR
mode=1;
tic
d1=fxyodrr(dn(:,:,:),flow,fhigh,dt,N,mode,verb,4);
toc
% figure;imagesc([d(:,:,9),dn(:,:,9),d1(:,:,9),dn(:,:,9)-d1(:,:,9)]);
% colormap(seis);
fprintf('SNR of RR is %g\n',yc_snr(d,d1,2));
%about 2.27s
%snr:4.49

%% DRR
mode=2;
tic
d2=fxyodrr(dn(:,:,:),flow,fhigh,dt,N,mode,verb,4);
toc
% figure;imagesc([d(:,:,9),dn(:,:,9),d2(:,:,9),dn(:,:,9)-d2(:,:,9)]);
% colormap(seis);
fprintf('SNR of DRR is %g\n',yc_snr(d,d2,2));
%about 2.13s %Note my version of optshrink is much faster (and even better) than the original one
%SNR: 10.13

%% ODRR
mode=3;
tic
d3=fxyodrr(dn(:,:,:),flow,fhigh,dt,N,mode,verb,4);
toc
% figure;imagesc([d(:,:,9),dn(:,:,9),d3(:,:,9),dn(:,:,9)-d3(:,:,9)]);
% colormap(seis);
fprintf('SNR of ODRR is %g\n',yc_snr(d,d3,2));
%about 2.45s
%SNR: 10.96

%% FK
tic
d4=yc_fkt(dn,'ps',10);
toc
fprintf('SNR of FK is %g\n',yc_snr(d,d4,2));
%about 0.015s
%SNR: 6.03


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



