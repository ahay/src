function Hyper(clean,noise,noisy,fxdecon,fxemd,fxemdpf)
%  Author      : Yangkang Chen
%                Texas Consortium of Computational Seismology
%		 Bureau of Economic Geology
%                Jackson School of Geosciences
%                The University of Texas at Austin
%         
%  Date        : Sep, 2013
%
%  Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
%    
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

% parameters definition
 flow=5;
 fhigh=120;
 dt=0.004;
 N=1;
 lf=15;
 mu=0.01;

% making events

   dt = 2./1000;
   tmax = 1.2;
   h = [-500:20:1000];
   tau = [0.1,.5,0.8];
   v = [1500,2400,2300];
   amp = [1., -1.,1];
   f0 = 20;
   snr = 2;
   L = 9;
   seed=2013;
d=hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L,seed); % noisy data
   snr = 200;
d0 = hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L,seed);  % clean data

n=d-d0;	%noise

% denoise
d2=fx_decon(d,dt,lf,mu,flow,fhigh);  			%using fx_decon
d3=fx_emd(d,flow,fhigh,dt,N);            		%using fx_emd
d4=fx_emdpf(d,flow,fhigh,dt,N, lf, mu);			%using fx_emdpf + 1 IMF

% from Matlab to Madagascar
rsf_create(clean,size(d0)');
rsf_write(d0,clean);

rsf_create(noise,size(n)');
rsf_write(n,noise);

rsf_create(noisy,size(d)');
rsf_write(d,noisy);

rsf_create(fxdecon,size(d2)');
rsf_write(d2,fxdecon);

rsf_create(fxemd,size(d3)');
rsf_write(d3,fxemd);

rsf_create(fxemdpf,size(d4)');
rsf_write(d4,fxemdpf);


















