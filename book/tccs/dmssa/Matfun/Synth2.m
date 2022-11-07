 function Synth2(clean,noisy,mssa,dmssa,K,N)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Apr, 2015

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
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

randn('state',201415);
shot0=shot;
n=0.5*randn(300,20,20);n=reshape(n,300,20*20);
n=yc_bp(n,0.002,1,5,120,130);
n=reshape(n,300,20,20);
shot=shot+n;%adding noise

%% Noise attenuation
d1=fxymssa(shot,0,120,0.004,K,0);	%MSSA
d2=fxydmssa(shot,0,120,0.004,K,N,0);	%DMSSA

%% 3D to 2D
shot0=reshape(shot0,300,20*20);
shot=reshape(shot,300,20*20);
d1=reshape(d1,300,20*20);
d2=reshape(d2,300,20*20);


%% from Matlab to Madagascar
rsf_create(clean,size(shot0)');
rsf_write(shot0,clean);

rsf_create(noisy,size(shot)');
rsf_write(shot,noisy);

rsf_create(mssa,size(d1)');
rsf_write(d1,mssa);

rsf_create(dmssa,size(d2)');
rsf_write(d2,dmssa);



