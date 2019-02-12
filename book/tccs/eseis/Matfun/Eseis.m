function Eseis(flag,infile1,n1,n2,N,flow,fhigh,dt,verb,outfile1,outfile2,infile2)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%
% Date        : Jan, 2015

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API

%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307

%% default values

if flag==1 % if flag=1: forward Eseis transform
    
    D=zeros(n1,n2);
    rsf_read(D,infile1);
    
    
    nf=2^nextpow2(n1);
    nf2=nf/2+1;
    
    % Transform into F-X domain
    DATA_FX=fft(D,nf,1);
    DATA_FX_RE=zeros(nf2,n2*N);
    DATA_FX_IM=zeros(nf2,n2*N);
    
    % First and last samples of the DFT.
    ilow  = floor(flow*dt*nf)+1;
    
    if ilow<1;
        ilow=1;
    end;
    
    ihigh = floor(fhigh*dt*nf)+1;
    
    if ihigh > nf2;
        ihigh=nf2;
    end
    
    % F-X domain emd
    for k=ilow:ihigh
        re=real(DATA_FX(k,:));
        im=imag(DATA_FX(k,:));
        imfre=emd(re,'MAXMODES',N-1);
        imfim=emd(im,'MAXMODES',N-1);
        [mr,nr]=size(imfre);
        [mi,ni]=size(imfim);
        if(mr<N || nr<n2) imfre=[imfre;zeros(N-mr,nr)]; imfre=[imfre zeros(N,n2-nr)]; end
        if(mi<N || ni<n2) imfim=[imfim;zeros(N-mi,ni)]; imfim=[imfim zeros(N,n2-ni)]; end

        DATA_FX_RE(k,:)=reshape(imfre.',1,N*n2);
        DATA_FX_IM(k,:)=reshape(imfim.',1,N*n2);
        
        if(mod(k,5)==0 && verb==1)
            fprintf( 'F %d is done!\n\n',k);
        end
    end
    
    
    
    
    %% from Matlab to Madagascar
    rsf_create(outfile1,size(DATA_FX_RE)');
    rsf_write(DATA_FX_RE,outfile1);
    
    rsf_create(outfile2,size(DATA_FX_IM)');
    rsf_write(DATA_FX_IM,outfile2);
    
else % if flag=0: Inverse Eseis transform
    
    nf=2^nextpow2(n1);
    nf2=nf/2+1;

    DATA_FX_RE=zeros(nf2,n2*N);
    DATA_FX_IM=zeros(nf2,n2*N);
    
    rsf_read(DATA_FX_RE,infile1);
    rsf_read(DATA_FX_IM,infile2);
   
    % Transform into F-X domain
    DATA_FX0=zeros(nf,n2);
    
    % First and last samples of the DFT.
    ilow  = floor(flow*dt*nf)+1;
    
    if ilow<1;
        ilow=1;
    end;
    
    ihigh = floor(fhigh*dt*nf)+1;
    
    if ihigh > nf2;
        ihigh=nf2;
    end
    
    for k=ilow:ihigh
        
        imfre=reshape(DATA_FX_RE(k,:),n2,N);
        imfim=reshape(DATA_FX_IM(k,:),n2,N);
        
        DATA_FX0(k,:)=sum(imfre,2)+i*sum(imfim,2);
        
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
    D1=D1(1:n1,:);
    
    rsf_create(outfile1,size(D1)');
    rsf_write(D1,outfile1);
    
end








