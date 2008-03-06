function BuildPareto(fname,xnname,rnname,xbnname,rbnname,xinname,rinname)
% Author      : G. Hennenfent
%               Seismic Laboratory for Imaging and Modeling
%               Department of Earth & Ocean Sciences
%               The University of British Columbia
%         
% Date        : March, 08

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2007 The University of British Columbia at Vancouver
%  Copyright (C) 2007 Gilles Hennenfent
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

jump = 6;
pts  = 1001;

load(fname)

xNorm = xNorm(:);
rNorm = rNorm(:);
gNorm = gNorm(:);

xbNorm = xNorm(1:jump:end);
rbNorm = rNorm(1:jump:end);
gbNorm = gNorm(1:jump:end);

xi = linspace(0,max(xNorm),pts);
yi = interp_qcmixed(xbNorm,rbNorm,-gbNorm./rbNorm,xi,0);
xi = xi(:);
yi = yi(:);

% unload xnorms
rsf_create(xnname,size(xNorm)');
rsf_write(xNorm,xnname);

% unload rnorms
rsf_create(rnname,size(rNorm)');
rsf_write(rNorm,rnname);

% unload xbnorms
rsf_create(xbnname,size(xbNorm)');
rsf_write(xbNorm,xbnname);

% unload rbnorms
rsf_create(rbnname,size(rbNorm)');
rsf_write(rbNorm,rbnname);

% unload xinorms
rsf_create(xinname,size(xi)');
rsf_write(xi,xinname);

% unload rinorms
rsf_create(rinname,size(yi)');
rsf_write(yi,rinname);
