function [d_interp] = spitz_fx_interpolation(d,dt,npf,pre1,pre2,flow,fhigh)
%SPITZ_FX_INTERPOLATION: Program for first order f-x interpolation using Spitz' method
%
%  [d_interp] = spitz_fx_interpolation(d,dt,npf,mu1,mu2,flow,fhigh)
%
%  IN   d:        data (d(nt,nh))
%       dt:       sampling interval in secs
%       npf:      prediction filter length
%       pre1:     percentage of pre-whitening for Pef estimation
%       pre2:     percentage of pre-whitening for estimation of missing samples
%       flow:     minimun temporal frequency to reconstruct (Hz)
%       fhigh:    maximun temporal frequency to reconstruct (Hz) 
%
%  OUT  d_iterp:  Interpolated data (d_interp(nt,2*nh-1))
%
%  Reference: Spitz, S., 1991, Seismic trace interpolation in the F-X 
%             domain: Geophysics, 56, 785-794.
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: Mostafa Naghizadeh (mnaghi@phys.ualberta.ca)
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

 [nt,nh] = size(d);

 nfft = 2^(nextpow2(nt));

 % Transform data to FX domain...

  % for estimating prediction filters  

    DF1 = fft(d,2*nfft,1);  

  % and for reconstruction

    DF2 = fft(d,nfft,1);   


 ilow = floor(flow*dt*nfft)+1;
 ihigh = floor(fhigh*dt*nfft)+1;

 INTDF = zeros(nfft, 2*nh-1);

 for ia = ilow:ihigh;

  % Select frequency slices from frequencies f and f/2

    x1 = DF1(ia,:); 
    x2 = DF2(ia,:); 
    
  % Estime the PFs needed for frequency sample ia;

    PF = prediction_filter(x1,npf,pre1);
    
  % Interpolate spatial data at freq. sample ia

    [y] = interpolate_freq(x2,PF,pre2);
    
  % Replace interpolated frequency slice into final matrix

    INTDF(ia,:) = y;
 
end

% Honor symmetries

 INTDF(nfft/2+2:nfft,:)=conj(flipud(INTDF(2:nfft/2,:)));

% Transform from f-x to t-x

 d_interp = real(ifft(INTDF,nfft,1));

 d_interp = d_interp(1:nt,:);

 return;


function [y]=interpolate_freq(x,PF,pre)
%INTERPOLATE_FREQ: Function used by spitz_fx_interpolation for reconstruction
%                  of traces
%
%  IN   x:   spatial data to be recostruction at a given freq.
%       PF:  prediction error filter
%       pre: prewhitening in percentage
%
%  OUT  y:   reconstructed vector of data at a particular frequency
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: Mostafa Naghizadeh (mnaghi@phys.ualberta.ca)
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



 np = length(PF);
 nx = length(x);
 ny = 2*nx-1;       % length of interpolated singal

% Forward step using PFs for interpolation

 TMPF1 = [fliplr(PF.'),-1];
 W1 = convmtx(TMPF1,ny-np);

% Bacward step

 TMPF2 = conj(fliplr(TMPF1));
 W2 = convmtx(TMPF2,ny-np);

% Forward and backward combined in an augmented system

 WT = [W1;W2];

% Separation of  Matrix into  known and unknown parts

 A = WT(:,2:2:ny);
 B = -1*WT(:,1:2:ny)*x.';

% Least squares solution for missing data, prewhitening
% is added to guarantee stability

 R= A'*A;
 g = A'*B;


 mu = (pre/100.)*trace(R)/(nx-1);

 y1 =  (R+mu*eye(nx-1))\g;

 y = zeros(1,ny);
 y(1:2:ny)=x;
 y(2:2:ny)=y1.';

 return;

function PF=prediction_filter(VEC,np,pre);
%PREDICTION_FILTER: Function used by spitz_fx_interpolation to compute
%                   prediction filters
%
%  IN   VEC:   spatial data at a given frequency
%       np:    length of pef
%       pre:   prewhitening in percentage
%
%  OUT  PF:    prediction filter
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: Mostafa Naghizadeh (mnaghi@phys.ualberta.ca)
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


 ns = length(VEC);

 for j=1:ns-np  

    C(j,:)=VEC(j+np:-1:j);

 end

 % Make kernel for pef estimation (forward and backward prediction included)

 A = [C(:,2:np+1);conj(fliplr(C(:,1:np)))];

 % RHS of system of equaitons

 B = [C(:,1);conj(C(:,np+1))];


 % Solve A.PF = B using least-squares with pre-whitening

 R = A'*A;
 g = A'*B;

 mu = (pre/100.)*trace(R)/np;

 PF =  (R+mu*eye(np))\g;

 return;


