function [d,h,t] = linear_events(dt,f0,tmax,h,tau,p,amp,snr,L,seed);
% LINEAR_EVENTS: A program to generate data containing linear events.
%
% [d] = linear_events(dt,f0,tmax,h,tau,p,amp,snr,L);
%
% IN   dt:        sampling interval in secs
%      f0:        central freq. of a Ricker wavelet in Hz
%      tmax:      maximun time of the simulation in secs
%      h:         vector of desire offsets in meters
%      tau,p,amp: vectors of intercept, ray parameter 
%                 and amplitude of each linear event
%                 (p is in sec/m and tau in secs)
%      snr:       signal to noise ratio (max amplitude in the clean
%                 signal/max amplitude of the band-pass noise)
%        L:       The random noise is average over L samples
%                 to simulate band-pass noise (L=1 means no averaging)
%	seed:	  seed for generating pseudo-random number
%
% OUT  d:         data that consist of a superposition of linear events
%  
% Example:        [d,h,t] = linear_events; imagesc(h,t,d);
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
  tmax = 1.;
  h = [0:10:10*(50-1)];
  tau = [0.1,0.2,0.3,0.6];
  p = [0.0004,-0.0001,0.0001,-0.0001];
  amp = [1.2,-1.,1.,1];
  f0 = 40;
  snr = 8;
  L = 5;
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

 delay = dt*(floor(nw/2)+1);

 for ifreq=1:nfft/2+1
  w = 2.*pi*(ifreq-1)/nfft/dt;
   for k=1:n_events
    Shift = exp(-i*w*(tau(k)+h*p(k)-delay));
   D(ifreq,:) = D(ifreq,:) + amp(k)* W(ifreq)*Shift;
  end
 end

% w-domain symmetries

 for ifreq=2:nfft/2
  D(nfft+2-ifreq,:) = conj(D(ifreq,:));
 end 

 d = ifft(D,[],1);
 d = real(d(1:nt,:));

% My definition of snr = (Max Amp of Clean Data)/(Max Amp of Noise)

 dmax  = max(max(abs(d)));
 op = hamming(L);
 Noise = conv2(randn(size(d)),op,'same');

 Noisemax = max(max(abs(Noise)));

 d = d + Noise*(dmax/Noisemax)/snr;

 if nargout>1;
  t = (0:1:nt-1)*dt; 
 end;

 return;
