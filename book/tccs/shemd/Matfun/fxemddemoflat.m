 function fxemddemoflat(clean,noise,noisy,denoised,noise1,fclean,fnoise,fnoisy,fdenoised,fnoise1)
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
[d] = linear_events(0.004,30,1.0, 1:100, [0.5],[0],1,20000,5,2013);  % dipping

d=d/max(max(d));

%rng(201313);
randn('state',201313);
n=0.1*randn(size(d)); %noise

D=d+n;


% parameters definition
 flow=5;
 fhigh=120;
 dt=0.004;
 N1=4;
 verb=1;
 Nf=20;

 %% Using fx_emd 
 d1= fx_emd(D,flow,fhigh,dt,N1,verb);  			
 n1=D-d1;
 %d12= fx_emd(D1,flow,fhigh,dt,2,verb);  			
 
% d13= fx_emd(D1,flow,fhigh,dt,3,verb);  
 
% d14= fx_emd(D1,flow,fhigh,dt,4,verb);
 
dfft=fft(d);
nfft=fft(n);
Dfft=fft(D);
d1fft=fft(d1);
n1fft=fft(n1);

fd1=real(dfft(Nf,:));           % frequency slice for clean data
fd2=real(nfft(Nf,:));           % frequency slice for noise
fd3=real(Dfft(Nf,:));           % frequency slice for noisy data
fd4=real(d1fft(Nf,:));          % frequency slice for denoised data
fd5=real(n1fft(Nf,:));          % frequency slice for denoised noise


 rsf_create(clean,size(d)');
 rsf_write(d,clean);

 rsf_create(noise,size(n)');
 rsf_write(n,noise);
 
 rsf_create(noisy,size(D)');
 rsf_write(D,noisy);

 rsf_create(denoised,size(d1)');
 rsf_write(d1,denoised);
 
 rsf_create(noise1,size(n1)');
 rsf_write(n1,noise1);
 
  rsf_create(fclean,size(fd1)');
 rsf_write(fd1,fclean);

 rsf_create(fnoise,size(fd2)');
 rsf_write(fd2,fnoise);
 
 rsf_create(fnoisy,size(fd3)');
 rsf_write(fd3,fnoisy);

 rsf_create(fdenoised,size(fd4)');
 rsf_write(fd4,fdenoised);
 
 rsf_create(fnoise1,size(fd5)');
 rsf_write(fd5,fnoise1);
 
