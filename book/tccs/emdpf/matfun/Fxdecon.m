function Fxdecon(flat,flat_fxdecon,dip,dip_fxdecon,complex,complex_fxdecon,complex_fxemdpf)
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
d1  =    linear_events(0.004,40,2,[0:10:10*79],1,0,1,2,2,seed);
%d10 =   linear_events(0.004,40,2,[0:10:10*79],1,0,1,200,2,seed);

d2=linear_events(0.004,40,2,[0:10:10*79],0.5,0.00126,1,2,2,seed);
%d20 = linear_events(0.004,40,2,[0:10:10*79],0.5,0.00126,1,200,2,seed);

d3 = linear_events(0.004,40,2,[0:10:10*79],[0.25,0.5, 0.85,1.25],[0,0.00126,0.0003,0.0006],[1,1,1,1],2,2,seed);
%d30 = linear_events(0.004,40,2,[0:10:10*79],[0.25,0.5, 0.85,1.25],[0,0.00126,0.0003,0.0006],[1,1,1,1],200,2,seed);

d11=fx_decon(d1,dt,lf,mu,flow,fhigh);
d22=fx_decon(d2,dt,lf,mu,flow,fhigh);
d33=fx_decon(d3,dt,lf,mu,flow,fhigh);
d44=fx_emdpf(d3,flow,fhigh,dt,N, lf,mu);

% from Matlab to Madagascar
%rsf_create(flat,size(d1)');
rsf_write(d1,flat);

%rsf_create(dip,size(d2)');
rsf_write(d2,dip);

%rsf_create(complex,size(d3)');
rsf_write(d3,complex);

%rsf_create(flat_fxdecon,size(d11)');
rsf_write(d11,flat_fxdecon);

%rsf_create(dip_fxdecon,size(d22)');
rsf_write(d22,dip_fxdecon);

%rsf_create(complex_fxdecon,size(d33)');
rsf_write(d33,complex_fxdecon);

%rsf_create(complex_fxemdpf,size(d44)');
rsf_write(d44,complex_fxemdpf);


