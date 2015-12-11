 function fxemddemodip(denoise,noise,fdenoise,fnoise)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Dec, 2013

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2013 The University of Texas at Austin
%  Copyright (C) 2013 Yangkang Chen
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
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% create synthetic data
[d] = linear_events(0.004,30,1.0, 1:100, [0.25],[0.005],1,20000,5,2013);  % dipping

d=d/max(max(d));

%rng(201313);
randn('state',201313);
n=0.1*randn(size(d)); %noise

D=d+n;


% parameters definition
 flow=5;
 fhigh=120;
 dt=0.004;
 N1=1;
 N2=2;
 N3=3;
 verb=1;
 Nf=20;

 %% Using fx_emd 
 d1= fx_mssa(D,flow,fhigh,dt,N1,verb);  			
 n1=D-d1;

 
dfft=fft(d1);
nfft=fft(n1);


fd1=real(dfft(Nf,:));           % frequency slice for denoised data
fd2=real(nfft(Nf,:));           % frequency slice for noise


 rsf_create(denoise,size(d1)');
 rsf_write(d1,denoise);

 rsf_create(noise,size(n1)');
 rsf_write(n1,noise);
 
 rsf_create(fdenoise,size(fd1)');
 rsf_write(fd1,fdenoise);
 
 rsf_create(fnoise,size(fd2)');
 rsf_write(fd2,fnoise);








 
