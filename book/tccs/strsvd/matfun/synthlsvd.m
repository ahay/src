 function synthlsvd(noisy,denoised,Nt,Nx,R,Twin,Xwin)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Apr, 2014

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
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% create synthetic data

nt=Nt;
nx=Nx;
param.r=R;
param.mode=1;
twin=Twin;
xwin=Xwin;

% allocate memory
data = zeros(nt,nx);

% from Madagascar to Matlab
rsf_read(data,noisy);

% processing in windows 
dlsvd=process_win(@localsvd,param,data,twin,xwin);

 
% from Matlab to Madagascar
rsf_create(denoised,size(dlsvd)');
rsf_write(dlsvd,denoised);

 
 
