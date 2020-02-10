 function Gentrace(t1,t2,dt,nt,f0,Q,trace)
% Author      : Yangkang Chen
%               Zhejiang University
%         
% Date        : Oct, 2018
%
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%
%  Copyright (C) 2018 Yangkang Chen
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
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% create synthetic trace with attenuation

t=0:dt:(nt-1)*dt;
nf=nt;
df=1.0/dt/(nf-1);
k1=t1/dt+1;
k2=t2/dt+1;

% f0 Hz Ricker wavelet
w1=yc_ricker(f0,dt);
nw1=size(w1,1);

% apply attenuation to the spectrum of Ricker wavelet
w2=ifft(fft([w1;zeros(nf-nw1,1)]).*exp(-[[0:fix(nf/2)]*df,[fix(nf/2)-1:-1:1]*df]'*pi*(t2-t1)/Q));
w2=w2(1:nw1);

% convolution
r1=zeros(nt,1);r1(k1)=10;
x1=conv(r1,w1,'same');

r2=zeros(nt,1);r2(k2)=10;
x2=conv(r2,w2,'same');

% normalize the data
x=x1+x2;x=x/max(x);

% from Matlab to Madagascar
rsf_create(trace,size(x)');
rsf_write(x,trace);
