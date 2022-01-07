function Synth(clean,noisy,obs,mssa,dmssa1,dmssa2,dmssa3)
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

%% without noise
dn=d;

%% decimate
[nt,nx,ny]=size(d);
ratio=0.5;
mask=genmask(reshape(d,nt,nx*ny),ratio,'c',201415);
mask=reshape(mask,nt,nx,ny);
d0=dn.*mask;

%% reconstruct
flow=0;fhigh=125;dt=0.004;N=3;Niter=10;mode=0;verb=1;
d1=fxymssa_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode);
% figure;imagesc([d(:,:,9),d0(:,:,9),d1(:,:,9),d(:,:,9)-d1(:,:,9)]);

%% simultaneous denoising and reconstruction
% adding noise
randn('state',201314);
var=0.2;
dn=d+var*randn(size(d));
d0=dn.*mask;

flow=0;fhigh=125;dt=0.004;N=3;Niter=10;mode=1;verb=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d1=fxymssa_recon(d0,mask,flow,fhigh,dt,N,Niter,eps,verb,mode,a);
% figure;imagesc([d(:,:,9),d0(:,:,9),d1(:,:,9),d(:,:,9)-d1(:,:,9)]);

% flow=0;fhigh=125;dt=0.004;N=3;Niter=10;mode=1;verb=1;NN=3;
% a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
% d2=fxydmssa_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,a);
% figure;imagesc([d(:,:,9),d0(:,:,9),d2(:,:,9),d(:,:,9)-d2(:,:,9)]);

% yc_snr(d,d0,2)
% yc_snr(d,d1,2)
% yc_snr(d,d2,2)
% SNR results
% -5.9853
% 1.4503
% 4.1965

%% from Matlab to Madagascar
rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);

rsf_create(obs,size(d0)');
rsf_write(d0,obs);

rsf_create(mssa,size(d1)');
rsf_write(d1,mssa);

flow=0;fhigh=125;dt=0.004;N=3;Niter=10;mode=1;verb=1;NN=1;
a=(Niter-(1:Niter))/(Niter-1); %linearly decreasing
d2=fxydmssa_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,a);

rsf_create(dmssa1,size(d2)');
rsf_write(d2,dmssa1);

NN=2;
d2=fxydmssa_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,a);

rsf_create(dmssa2,size(d2)');
rsf_write(d2,dmssa2);

NN=3;
d2=fxydmssa_recon(d0,mask,flow,fhigh,dt,N,NN,Niter,eps,verb,mode,a);

rsf_create(dmssa3,size(d2)');
rsf_write(d2,dmssa3);



