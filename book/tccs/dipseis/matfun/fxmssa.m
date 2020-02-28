function [ D1 ] = fxmssa(D,flow,fhigh,dt,N,verb)
%  FX_MSSA: F-X domain multichannel singular spectrum analysis (MSSA)
%
%  IN   D:   	intput data 
%       flow:   processing frequency range (lower)
%       fhigh:  processing frequency range (higher)
%       dt:     temporal sampling interval
%       N:      number of singular value to be preserved
%       verb:   verbosity flag (default: 0)
%      
%  OUT  D1:  	output data
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
%  Reference:   Simultaneous seismic data denoising and reconstruction via multichannel 
%               singular spectrum analysis, Geophysics, 2011, 76, V25-V32
%

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

lx=floor(trace/2)+1;
lxx=trace-lx+1;

% F-X domain emd
for k=ilow:ihigh
   
    r=hankel(DATA_FX(k,1:lx),[DATA_FX(k,lx:trace)]);

    [U,E,V]=svd(r);
    
%    D1=D; % the size of D1 is round((trace+1)/2)*trace+1-round((trace+1)/2)
%    for j=N+1:trace+1-round((trace+1)/2)
%        D1(j,j)=0;
%    end
    r1=U(:,1:N)*E(1:N,1:N)*V(:,1:N)';   
%     DATA_FX0(k,:)=[r1(:,1).',r1(round((trace+1)/2),2:trace+1-round((trace+1)/2))];
    DATA_FX0(k,:)=ave_antid(r1);
     
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

function [dout] =ave_antid(din);
% averaging along antidiagonals

   [n1,n2]=size(din);
   nout=n1+n2-1;
   dout=zeros(nout,1);
   for i=1:nout	
       if i<n1
          for id=1:i
	  	    dout(i)=dout(i) + din(i-(id-1),1+(id-1))/i; 
	  	end
	  else
          for id=1:nout+1-i
	  	    dout(i)=dout(i) + din(n1-(id-1),1+(i-n1)+(id-1))/(nout+1-i); 
	  	end	  
	  end
   end
return
