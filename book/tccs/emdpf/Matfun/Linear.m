function Linear(clean,noise,noisy,fxdecon,fxemd,fxemdpf1,fxemdpf2,fxemdpf3)
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
 seed=2013;



% making events
[d1,h,t] = linear_events(0.004,40,2,1:40,0.5,0,2,2,1,seed);
[d2,h,t] = linear_events(0.004,40,2,1:40,0.5,0.025,2,2,1,seed);
[d3,h,t] = linear_events(0.004,40,2,1:40,1.5,0,2,2,1,seed);
 
d=[d1,d2,d3];     % noisy data
 
[d10,h,t] = linear_events(0.004,40,2,1:40,0.5,0,2,200,1,seed);
[d20,h,t] = linear_events(0.004,40,2,1:40,0.5,0.025,2,200,1,seed);
[d30,h,t] = linear_events(0.004,40,2,1:40,1.5,0,2,200,1,seed);
 
d0=[d10,d20,d30]; % clean data

n=d-d0;		 % noise

% denoise
d2=fx_decon(d,dt,lf,mu,flow,fhigh);  			%using fx_decon
d3=fx_emd(d,flow,fhigh,dt,N);            		%using fx_emd
d4=fx_emdpf(d,flow,fhigh,dt,N, lf, mu);			%using fx_emdpf + 1 IMF
d5=fx_emdpf(d,flow,fhigh,dt,N+1, lf, mu);		%using fx_emdpf + 2 IMF
d6=fx_emdpf(d,flow,fhigh,dt,N+2, lf, mu);		%using fx_emdpf + 3 IMF

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

rsf_create(fxemdpf1,size(d4)');
rsf_write(d4,fxemdpf1);
rsf_create(fxemdpf2,size(d5)');
rsf_write(d5,fxemdpf2);
rsf_create(fxemdpf3,size(d6)');
rsf_write(d6,fxemdpf3);

















