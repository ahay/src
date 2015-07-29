function [d,h,t] = hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L,seed);
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
%	seed:	  seed for generating pseudo-random number
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
%  Copyright (C) 2013, Texas Consortium for computational seismology
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi, Yangkang Chen
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
  seed=2013;
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

rng(seed);
noisetemp=randn(size(d));
 Noise = conv2(noisetemp,op,'same');

 Noisemax  = max(max(Noise));

 d = d + (dmax/Noisemax)*Noise/snr;

 if nargout>1;
  t = (0:1:nt-1)*dt;
 end;

 return;
