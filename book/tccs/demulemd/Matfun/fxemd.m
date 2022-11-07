function [ D1 ] = fxemd(D,flow,fhigh,dt,N,verb)
%FXEMD: F-X domain empirical mode decomposition along the spatial dimension
%  IN   D:   	intput data 
%       flow:  processing frequency range (lower)
%       fhigh: processing frequency range (higher)
%       dt:    temporal sampling interval
%       N: 	number of IMF to be removed
%       verb:   verbosity flag (default: 0)
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
%  REFERENCES:
%  Chen, Y. and J. Ma, 2014, Random noise attenuation by f-x empirical mode decomposition predictive filtering, Geophysics, 79, V81-V91.
%  Chen, Y., C. Zhou, J. Yuan, and Z. Jin, 2014, Application of empirical mode decomposition to random noise attenuation of seismic data, Journal of Seismic Exploration, 23, 481-495.
%  Chen, Y., S. Gan, T. Liu, J. Yuan, Y Zhang, and Z. Jin, 2015, Random noise attenuation by a selective hybrid approach using f-x empirical mode decomposition, Journal of Geophysics and Engineering, 12, 12-25.
%  Chen, Y., G. Zhang, S. Gan, and C. Zhang, 2015, Enhancing seismic reflections using empirical mode decomposition in the flattened domain, Journal of Applied Geophysics, 119, 99-105.
%  Gan, S., S. Wang, Y. Chen, J. Chen, W. Zhong, and C. Zhang, 2016, Improved random noise attenuation using f − x empirical mode decomposition and local similarity, Applied Geophysics, 13, 2016, 127-134.
%  Chen, Y., 2016, Dip-separated structural filtering using seislet thresholding and adaptive empirical mode decomposition based dip filter, Geophysical Journal International, 206, 457-469.
%  Chen, W., J. Xie, S. Zu, S. Gan, and Y. Chen, 2017, Multiple reflections noise attenuation using adaptive randomized-order empirical mode decomposition, IEEE Geoscience and Remote Sensing Letters, 14, 18-22.
%  Chen, Y., Y. Zhou, W. Chen, S. Zu, W. Huang, and D. Zhang, 2017, Empirical low rank decomposition for seismic noise attenuation, IEEE Transactions on Geoscience and Remote Sensing, 2017, 55, 4696-4711.
%  Chen, Y. and S. Fomel, 2018, EMD-seislet transform, Geophysics, 83, A27-A32.

if nargin==0
 error('Input data must be provided!');
end

if nargin==1
 flow=1;
 fhigh=124;
 dt=0.004;
 N=1;
 verb=0;
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

% F-X domain emd
for k=ilow:ihigh
    re=real(DATA_FX(k,:));
    im=imag(DATA_FX(k,:));
    imfre=emd(re,'MAXMODES',N);
    imfim=emd(im,'MAXMODES',N);
    [mr,nr]=size(imfre);
    [mi,ni]=size(imfim);
    if(mr<N || nr<trace) imfre=[imfre;zeros(N-mr,nr)]; imfre=[imfre zeros(N,trace-nr)]; end
    if(mi<N || ni<trace) imfim=[imfim;zeros(N-mi,ni)]; imfim=[imfim zeros(N,trace-ni)]; end
           
    if(N==1)
        DATA_FX0(k,:)=re-imfre(1:N,:)+i*(im-imfim(1:N,:));
    else
        DATA_FX0(k,:)=re-sum(imfre(1:N,:))+i*(im-sum(imfim(1:N,:)));
    end 
    if(mod(k,5)==0 && verb==1)
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
