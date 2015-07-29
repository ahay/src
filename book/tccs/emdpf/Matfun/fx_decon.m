function [DATA_f] = fx_decon(DATA,dt,lf,mu,flow,fhigh);
%FX_DECON: SNR enhancement using fx-deconvolution.
%
%  [DATA_f] = fx_decon(DATA,dt,lf,mu,flow,fhigh);
% 
%  IN   DATA:   the data matrix, columns are traces
%       dt:     sampling interval in sec
%       lf:     lenght of operator (lenght of the filter)
%       mu:     pre-whitening 
%       flow:   min  freq. in the data in Hz
%       fhigh:  max  freq. in the data in Hz
% 
%  OUT  DATA_f: filtered data 
%
% 
%  Reference: Canales, 1984, Random noise reduction, 54.th. Ann. Internat. 
%             Mtg., Soc. Expl. Geophys., Expanded Abstracts, pp. 525-527
%
%  Note: Canales method is modified to use non-Toeplitz system of equations
%        with backward and foward prediction filters
%       
%  Example: see fx_decon_demo.m
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
%


 [nt,ntraces] = size(DATA);
 nf = 2^nextpow2(nt);
 
 DATA_FX_f = zeros(nf,ntraces);
 DATA_FX_b = zeros(nf,ntraces);

% First and last samples of the DFT.

 ilow  = floor(flow*dt*nf)+1; 

  if ilow<1; 
   ilow=1; 
  end;

 ihigh = floor(fhigh*dt*nf)+1;

  if ihigh > floor(nf/2)+1; 
   ihigh=floor(nf/2)+1; 
  end

% Transform to FX

 DATA_FX = fft(DATA,nf,1);

 for k = ilow:ihigh;
  aux_in  = DATA_FX(k,:)';
  [aux_out_f,aux_out_b] = ar_modeling(aux_in,lf,mu);
  DATA_FX_f(k,:) = aux_out_f';
  DATA_FX_b(k,:) = aux_out_b';
 end;

% Honor symmetries

 for k=nf/2+2:nf
  DATA_FX_f(k,:) = conj(DATA_FX_f(nf-k+2,:));
  DATA_FX_b(k,:) = conj(DATA_FX_b(nf-k+2,:));
 end

% Back to TX (the output) 

 DATA_f = real(ifft(DATA_FX_f,[],1));
 DATA_f = DATA_f(1:nt,:);

 DATA_b = real(ifft(DATA_FX_b,[],1));
 DATA_b = DATA_b(1:nt,:);

% Average predictions (forward and backward)

 DATA_f = (DATA_f + DATA_b);
 DATA_f(:,lf+1:ntraces-lf)= DATA_f(:,lf+1:ntraces-lf)/2;
 
return

function [yf,yb] = ar_modeling(x,lf,mu);
%AR_MODELING: autoregressive modeling of 1D spatial data
%
%  IN    x:   data 
%        lf:  length of the operator
%        mu:  pre-whitening in %
%      
%  OUT   yf:  prediction of the data using forward AR modeling
%        yb:  prediction of the data using backward AR modeling
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



   nx = length(x);

% backward ar-modeling

   y  = x(1:nx-lf,1);
   C  = x(2:nx-lf+1,1);
   R  = x(nx-lf+1:nx,1);
   M = hankel(C,R);

   B = M'*M;  beta = B(1,1)*mu/100;
   ab = (B + beta*eye(lf))\M'*y;
   temp = M*ab;
   temp = [temp;zeros(lf,1)];
   yb = temp;


   y  = x(lf+1:nx,1);
   C  = x(lf:nx-1,1);
   R = flipud(x(1:lf,1));
   M = toeplitz(C,R);


   B = M'*M;  beta = B(1,1)*mu/100;

   af = (B + beta*eye(lf))\M'*y;
   temp = M*af;
   temp = [zeros(lf,1);temp];
   yf = temp;

return
