function [ D1 ] = fx_emdpf(D,flow,fhigh,dt,N, lf, mu)
%FXEMDPF: F-X domain empirical mode decomposition predictive filtering 
%  IN    D:   	intput data 
%        flow:  processing frequency range (lower)
%        fhigh: processing frequency range (higher)
%	 dt:    temporal sampling interval
%	 N: 	number of IMF to be extracted for predictive filtering
%	 lf:    predictive length
%	 mu:    regularization term
%      
%  OUT   D1:  	output data
%
%  Copyright (C) 2013 The University of Texas at Austin
%  Copyright (C) 2013 Yangkang Chen
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
if nargin==0
 error('Input data must be provided!');
end

if nargin==1
 flow=0;
 fhigh=125;
 dt=0.004;
 N=1;
 lf=15;
 mu=0.01;
end;

[sample,trace ]=size(D);
D1=zeros(sample,trace);

nf=2^nextpow2(sample);
nk=2^nextpow2(trace);

% Transform into F-X domain
DATA_FX=fft(D,nf,1);
DATA_FX0=zeros(nf,trace);

% First and last samples of the DFT.
ilow  = floor(flow*dt*nf)+1;

if ilow<1;
    ilow=1;
end;

ihigh = floor(fhigh*dt*nf)+1;

if ihigh > floor(nf/2)+1;
    ihigh=floor(nf/2)+1;
end

% F-X domain emdpf
for k=ilow:ihigh
    re=real(DATA_FX(k,:));
    im=imag(DATA_FX(k,:));
    imfre=emd(re,'MAXMODES',N);
    imfim=emd(im,'MAXMODES',N);
    if(N==1)
        renew=re-imfre(1:N,:);
	imnew=im-imfim(1:N,:);
	DATA_FX0(k,:)=renew+i*imnew;
	[yf,yb] = ar_modeling((imfre(1:N,:)+i*imfim(1:N,:))',lf,mu);
    else
        renew=re-sum(imfre(1:N,:));
	imnew=im-sum(imfim(1:N,:));
	DATA_FX0(k,:)=renew+i*imnew;
	[yf,yb] = ar_modeling((sum(imfre(1:N,:))+i*sum(imfim(1:N,:)))',lf,mu);
    end 
    y_ave = (yf' + yb');
    y_ave(:,lf+1:trace-lf)= y_ave(:,lf+1:trace-lf)/2;
    DATA_FX0(k,:)=DATA_FX0(k,:)+y_ave;
    if(mod(k,5)==0)
        fprintf( 'F %d is done!\n\n',k);
    end
end

% Honor symmetries
for k=nf/2+2:nf
    DATA_FX0(k,:) = conj(DATA_FX0(nf-k+2,:));
end

% Back to TX (the output)
D1=real(ifft(DATA_FX0,[],1));
D1=D1(1:sample,:);

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
