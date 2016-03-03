function Syn_ddtf(din,dout,n1,n2,lambda,niter,lvl,htype,thr,thrtype,sh)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : May, 2014

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
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
if nargin==9
	thrtype='e';
	sh='s';
end;

if nargin==10
	sh='s';
end;

un=zeros(n1,n2); % noisy input
rsf_read(un,din);

%% reference parameters
%    lambda=0.5;
%    niter=30;
%    lvl=2;
%    htype='spline';
%    thr=0.02; or 20


u1=ddtf_denoise2d(un, lambda, niter, lvl, htype, thr, thrtype, sh);

% from Matlab to Madagascar
rsf_create(dout,size(u1)');
rsf_write(u1,dout);




















