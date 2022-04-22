function [ D1 ] = fx_emd(D,flow,fhigh,dt,N)
%FXEMD: F-X domain empirical mode decomposition along the spatial dimension
%  IN    D:   	intput data 
%        flow:  processing frequency range (lower)
%        fhigh: processing frequency range (higher)
%	 dt:    temporal sampling interval
%	 N: 	number of IMF to be removed
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
    if(N==1)
        DATA_FX0(k,:)=re-imfre(1:N,:)+i*(im-imfim(1:N,:));
    else
        DATA_FX0(k,:)=re-sum(imfre(1:N,:))+i*(im-sum(imfim(1:N,:)));
    end 
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
