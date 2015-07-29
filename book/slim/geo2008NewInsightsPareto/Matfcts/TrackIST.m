function TrackIST(Aname,x0name,xnname,rnname,Iters,lambda)
% Author      : G. Hennenfent
%               Seismic Laboratory for Imaging and Modeling
%               Department of Earth & Ocean Sciences
%               The University of British Columbia
%         
% Date        : November, 07

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
    
dims = rsf_dim(Aname);
n = dims(1);
N = dims(2);

% load matrix
A = zeros(n,N);
rsf_read(A,Aname);

% load sparse vector
x0 = zeros(N,1);
rsf_read(x0,x0name);

% compute data
y = A*x0;

% compute path of ISTc for BP
[x, xnorms, rnorms] = SolveIST(A, y, Iters, lambda, [], 1);

% unload xnorms
rsf_create(xnname,[length(xnorms);1]);
rsf_write(xnorms,xnname);

% unload rnorms
rsf_create(rnname,[length(rnorms);1]);
rsf_write(rnorms,rnname);
