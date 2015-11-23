function [ D1 ] = fx_emdmssa(D,flow,fhigh,dt,N1, N2, verb)
%FXEMDPF: F-X domain empirical mode decomposition combined with multichannel singular spectrum analysis (MSSA)
%
%  IN   D:   	intput data 
%       flow:  processing frequency range (lower)
%       fhigh: processing frequency range (higher)
%       dt:    temporal sampling interval
%       N1: 	number of IMF to be extracted for MSSA
%       N2:     number of singular value to be preserved 
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
 N1 =1;
 N2 =1;
 verb=1;
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
    imfre=emd(re,'MAXMODES',N1);
    imfim=emd(im,'MAXMODES',N1);
    [mr,nr]=size(imfre);
    [mi,ni]=size(imfim);
    if(mr<N1) imfre=[imfre;zeros(N1-mr,nr)]; end
    if(mi<N1) imfim=[imfim;zeros(N1-mi,ni)]; end
    if(N1==1)
        renew=re-imfre(1:N1,:);
        imnew=im-imfim(1:N1,:);
        DATA_FX0(k,:)=renew+i*imnew;
        
        d_mssa= imfre(1:N1,:)+i*imfim(1:N1,:);
        r=hankel(d_mssa(1:round((trace+1)/2)),[d_mssa(round((trace+1)/2):trace)]);
        [U,D,V]=svd(r);
        D1=D; % the size of D1 is round((trace+1)/2)*trace+1-round((trace+1)/2)
        for j=N2+1:trace+1-round((trace+1)/2)
            D1(j,j)=0;
        end
        r1=U*D1*V';    
        d_mssa=[r1(:,1).',r1(round((trace+1)/2),2:trace+1-round((trace+1)/2))];
    
    else
        renew=re-sum(imfre(1:N1,:));
        imnew=im-sum(imfim(1:N1,:));
        DATA_FX0(k,:)=renew+i*imnew;
  
        d_mssa= sum(imfre(1:N1,:))+i*sum(imfim(1:N1,:));
        r=hankel(d_mssa(1:round((trace+1)/2)),[d_mssa(round((trace+1)/2):trace)]);
        [U,D,V]=svd(r);
        D1=D; % the size of D1 is round((trace+1)/2)*trace+1-round((trace+1)/2)
        for j=N2+1:trace+1-round((trace+1)/2)
            D1(j,j)=0;
        end
        r1=U*D1*V';    
        d_mssa=[r1(:,1).',r1(round((trace+1)/2),2:trace+1-round((trace+1)/2))];
        
    end 
       
    DATA_FX0(k,:)=DATA_FX0(k,:)+d_mssa;
    
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


