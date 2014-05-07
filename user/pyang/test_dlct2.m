function test_dlct2

% Reference: Alkishriwo, Osama A., and Luis F. Chaparro. "A discrete 
%  linear chirp transform (DLCT) for data compression." Information 
%  Science, Signal Processing and their Applications (ISSPA), 2012 
%  11th International Conference on. IEEE, 2012.

clc,clear, close all


dt = 2./1000;
tmax = 1.2;
h = [-500:18:1200];
tau = [0.1,.35,0.8];
v = [1500,2400,2300];
amp = [1., -1.,1];
f0 = 20;
snr = 2;
Ln = 9;
d = hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,Ln);
d_free=hyperbolic_events(dt,f0,tmax,h,tau,v,amp,60,Ln);
[nt,nx]=size(d);


L=100;
C=0.01;
pclip=0.5;
% -------------decomposition-----------------------
u0=zeros(nt,L,nx);
u1=zeros(nt,L,nx);
d_new=zeros(nt,nx);
figure(1),clf
for ix=1:nx
    slice=dlct(d(:,ix), C, L,nt);
    u0(:,:,ix)=fft(slice,[],2);   
    subplot(211)
    imagesc(abs(slice))
    subplot(212)
    imagesc(abs(fftshift(u0(:,:,ix))))
    drawnow,%pause
end

%------------------- thresholding ---------------------
[b,ind]=sort(abs(u0));
thr=abs(u0(ind(floor(pclip*nt*L*nx))));
u1=thresh(u0,thr);

%------------------- reconstruction ---------------------
for ix=1:nx
    slice=ifft(u1(:,:,ix),[],2);
    d_new(:,ix)=idlct(slice,C,L,nt);  
end

figure(2),clf
subplot(211)
imagesc(d)
subplot(212)
imagesc(d_new)
snr1=mysnr(d,d_free)
snr2=mysnr(d_new, d_free)

function y=mysnr(x, xref)
err=x-xref;
ss=sum(abs(x(:)).^2);
nn=sum(abs(err(:)).^2);
y=10*log10(ss/nn);


function y=thresh(x,thr)
y=x.*(abs(x)>thr);




function y=dlct(x,C, L,N)
%  forward discrete linear chirp transform
n=1:N;
h=zeros(N,1);
y=zeros(N,L);
for l=-L/2:L/2-1
    h=x.*exp(-j*2*pi*C*l*n'.^2/N);
    y(:,l+L/2+1)=fft(h);
end

function x=idlct(y,C,L,N)
%  forward discrete linear chirp transform
n=1:N;
g=zeros(N,L);
for l=-L/2:L/2-1
    g(:,l+L/2+1)=ifft(y(:,l+L/2+1));
end
x=sum(g.*exp(j*2*pi*C*n'.^2*[-L/2:L/2-1]/N),2)/L;
x=real(x);


function [d,h,t] = hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L);
%HYPERBOLIC_EVENTS: A program to generate data containing hyperbolas.
%
%  [d] = hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L);
%
%  IN   dt:        sampling interval in secs
%       f0:        central freq. of a Ricker wavelet in Hz
%       tmax:      maximun time of the simulation in secs
%       h:         vector of offsets in meters
%       tau,v,amp: vectors of intercept, rms velocities
%                  and amplitude of each linear event
%                  (v is in m/s and tau in secs)
%       snr:       signal to noise ratio (max amplitude of the clean
%                  signal/max amplitude of the noise)
%       L:         The random noise is average over L samples
%                  to simulate band-pass noise (L=1 means no averaging)
%
%  OUT  d:         Data that consist of a superposition of reflections
%                  with hyerbolic  moveout (no avo)
%       t,h:       time and offset axes 
%
%  Example with default parameters:
%
%    [d,h,t] = hyperbolic_events; imagesc(h,t,d);
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%


 if nargin == 0
  dt = 2./1000;
  tmax = 1.2;
  h = [-500:20:1000];
  tau = [0.1,.5,0.8];
  v = [1500,2400,2300];
  amp = [1., -1.,1];
  f0 = 20;
  snr = 2;
  L = 9;
 end;
 
 nt = floor(tmax/dt)+1;
 nfft = 4*(2^nextpow2(nt));
 n_events = length(tau);
 nh = length(h);
 wavelet = ricker(f0,dt); 
 nw = length(wavelet);
 W = fft(wavelet,nfft);
 D = zeros(nfft,nh);
 i = sqrt(-1);

% Important: the following lines is to have the maximum of the Ricker
% wavelet at the right intercept time

 delay = dt*(floor(nw/2)+1);

 for ifreq=1:nfft/2+1
  w = 2.*pi*(ifreq-1)/nfft/dt;
   for k=1:n_events
    Shift = exp(-i*w*(  sqrt(tau(k)^2 + (h/v(k)).^2) - delay));
   D(ifreq,:) = D(ifreq,:) +amp(k)* W(ifreq)*Shift;
  end
 end

% Apply w-domain symmetries

 for ifreq=2:nfft/2
  D(nfft+2-ifreq,:) = conj(D(ifreq,:));
 end 

 d = ifft(D,[],1);
 d = real(d(1:nt,:));

 dmax  = max(max(d));
 op = hamming(L);
 ops = sum(sum(op));
 op = op/ops;
 Noise = conv2(randn(size(d)),op,'same');

 Noisemax  = max(max(Noise));

 d = d + (dmax/Noisemax)*Noise/snr;

 if nargout>1;
  t = (0:1:nt-1)*dt;
 end;

 return;
 
 
function [w,tw] = ricker(f,dt)
%RICKER: Ricker wavelet of central frequency f.
%
%  [w,tw] = ricker(f,dt);
%
%  IN   f : central freq. in Hz (f <<1/(2dt) )
%       dt: sampling interval in sec  
%
%  OUT  w:  the Ricker wavelet
%       tw: axis
%
%  Example
%
%    [w,tw] = ricker(10,0.004);
%    plot(tw,w);
%

%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%


 nw=2.2/f/dt;
 nw=2*floor(nw/2)+1;
 nc=floor(nw/2);
 w = zeros(nw,1);

 k=[1:1:nw]';

 alpha = (nc-k+1).*f*dt*pi;
 beta=alpha.^2;
 w = (1.-beta.*2).*exp(-beta);

  if nargout>1;
    tw = -(nc+1-[1:1:nw])*dt;
  end







