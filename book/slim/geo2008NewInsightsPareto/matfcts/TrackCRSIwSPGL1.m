function TrackCRSIwSPGL1(dname,iname,fname,xnname,rnname,maxiters)
% Author      : G. Hennenfent
%               Seismic Laboratory for Imaging and Modeling
%               Department of Earth & Ocean Sciences
%               The University of British Columbia
%         
% Date        : November, 07

% Requirements: SPGL1 (http://www.cs.ubc.ca/labs/scl/index.php/Main/Spgl1)
%               RSF (http://rsf.sourceforge.net/) with Matlab API
%               Sparco (http://www.cs.ubc.ca/labs/scl/sparco/)
    
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

P = generateProblem(901);
A = P.A;
y = P.b;
f = P.signal;

% set SPGL1 parameters
opts  = spgSetParms('iterations', maxiters, ...
                    'verbosity' ,        2, ...
                    'bpTol'     ,     1e-5, ...
                    'optTol'    ,     1e-5, ...
                    'decTol'    ,     1e-3, ...
                    'subspaceMin',    0     ...
                    );
sigma = 0;
tau   = 0;

% compute Pareto curve
[x, r, g, info] = spgl1(A, y, tau, sigma, [], opts);

interp = P.reconstruct(x);
data = reshape(P.op.Restrict(y,2),P.signalSize); 

% unload interp
rsf_create(iname,size(interp)');
rsf_write(interp-mean(f(:)),iname);

% unload model
rsf_create(fname,size(f)');
rsf_write(f-mean(f(:)),fname);

% unload data
rsf_create(dname,size(data)');
rsf_write(data-mean(f(:)),dname);

% unload xnorms
rsf_create(xnname,[length(info.xNorm1);1]);
rsf_write(info.xNorm1,xnname);

% unload rnorms
rsf_create(rnname,[length(info.rNorm2);1]);
rsf_write(info.rNorm2,rnname);
