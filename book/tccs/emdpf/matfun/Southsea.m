function Southsea(data,fxdecon,fxemd,fxemdpf)
%  Author      : Yangkang Chen
%                Texas Consortium of Computational Seismology
%		 Bureau of Economic Geology
%                Jackson School of Geosciences
%                The University of Texas at Austin
%         
%  Date        : Sep, 2013

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

% parameters definition
 flow=5;
 fhigh=245;
 dt=0.002;
 N=3;
 lf=10;
 mu=0.01;
 twin=500;
 xwin=500;
 mode1=1;
 mode2=2;
 mode3=3;

% allocate memory
d = zeros(2500,2000);

% from Madagascar to Matlab
rsf_read(d,data);

% denoise
d1= denoise_win(d,flow,fhigh,mode1,dt,N,lf,mu,twin,xwin);  			%using fx_decon
d2= denoise_win(d,flow,fhigh,mode2,dt,N,lf,mu,twin,xwin);  			%using fx_emd
d3= denoise_win(d,flow,fhigh,mode3,dt,N,lf,mu,twin,xwin);  			%using fx_emdpf


% from Matlab to Madagascar
rsf_create(data,size(d)');
rsf_write(d,data);

rsf_create(fxdecon,size(d1)');
rsf_write(d1,fxdecon);

rsf_create(fxemd,size(d2)');
rsf_write(d2,fxemd);

rsf_create(fxemdpf,size(d3)');
rsf_write(d3,fxemdpf);


















