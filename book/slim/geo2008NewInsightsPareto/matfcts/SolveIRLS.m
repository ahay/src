function [x, xnorms, rnorms] = SolveIRLS(A, y, Iters, InnerIters, a, sigma, damp, x, fullPath)
% SolveIRLS: Iterative Reweighted least-squares
%
%---------------------------------------------------------------
% Solve
%
%      min 1/2||y-Ax||_2^2 + damp*||Wx||_2^2
%
%---------------------------------------------------------------
%
% INPUTS
% ======
% A             is an n-by-N matrix, explicit or an operator.
%               If A is a function, then it must have the signature
% 
%               y = A(x,mode)   if mode == 1 then y = A x  (y is n-by-1);
%                               if mode == 2 then y = A'x  (y is N-by-1).
% y             is an m-vector.
% Iters         Total number of iterations
% InnerIters    Total number of iterations for LSQR
% a             is a tuning parameter for sparsity.
%               a=1/2 -> l1-2 (quadratic approximation of l1-norm)
%               a=1   -> Cauchy's norm
%               a=2   -> geman-McClure (sparsest)
% sigma         stabilization parameter (small -> sparser for fixed a)
% x             is an N-vector estimate of the solution (possibly all
%               zeros). If x = [], then ISTc determines N via
%               N = length( A'y ) and sets  x = zeros(N,1).
% fullPath      is 1 to track the path, 0 otherwise.
%
% OUTPUTS
% =======
% x             is a solution of the problem
% xnorms        is a MaxIters-vector of the l1-norm of the estimates
%               of the solution at each iteration.
% rnorms        is a MaxIters-vector of the l2-norm of the residual 
%               at each iteration.

% Author      : G. Hennenfent
%               Seismic Laboratory for Imaging and Modeling
%               Department of Earth & Ocean Sciences
%               The University of British Columbia
%         
% Date        : November, 07

%  Copyright (C) 2007 The University of British Columbia at Vancouver
%  Copyright (C) 2007 Gilles Hennenfent and Ewout van den Berg
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

%----------------------------------------------------------------------
% Check arguments. 
%----------------------------------------------------------------------
if nargin < 9
    fullPath = 0;
end

if nargin < 8
    x = [];
end

if nargin < 7
    damp = 1;
end

if nargin < 6
    sigma = [];
end

if nargin < 5
    a = .5;
end

if nargin < 4
    InnerIters = 100;
end

if nargin < 3
    Iters = 500;
end

if nargin < 2
    error('At least two arguments are required');
end

OuterIters = ceil(Iters/InnerIters);
iter = 0;

xnorms = [];
rnorms = [];

%----------------------------------------------------------------------
% Determine initial x and vector length N
%----------------------------------------------------------------------
explicit = isnumeric(A);
object   = isobject(A);

projres = Aprod(y,2);
N = length(projres);
M = length(y);

if isempty(x)
    x = zeros(N,1);
end

if isempty(sigma)
    sigma = sqrt(tsprctile(projres.^2,10));
end

rnorms(1) = norm(y-Aprod(x,1),2);
xnorms(1) = norm(x,1);

W = ones(N,1);
r = y;
z = zeros(N,1);
rd = [r; -damp*z];

%----------------------------------------------------------------------
% MAIN LOOP
%----------------------------------------------------------------------
h = waitbar(0,'SolveIRLS - Please wait...');
for i = 1:OuterIters
    % FINE!!!
    z = lsqr(prod(size(y)),N,@Mprod,y,damp,[],[],[],InnerIters,2);
    x = W.*z;
    r  = y - Aprod(x,1);
    
    % EXPERIMENTAL
% $$$     dz = lsqr(M+N,N,@MIprod,rd,0,[],[],[],InnerIters,2);
% $$$     z  = z + dz;
% $$$     x  = W.*z;
% $$$     r  = y - Aprod(x,1);
% $$$     rd = [r; -damp*z];

    if fullPath
         rnorms(i+1) = norm(r,2);
         xnorms(i+1) = norm(x,1);
    end
    iter = iter+1;
    waitbar(iter/OuterIters,h)
    if iter == OuterIters
        close(h)
        break
    else
        W = (x.*conj(x) + sigma^2).^(a/2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS.  These share some vars with workspace above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function z = Aprod(x,mode)
   if mode == 1
      if   explicit || object
          z = A*x;
      else 
          z = A(x,1);
      end
   elseif mode == 2
      if   explicit || object
          z = A'*x;
      else
          z = A(x,2);
      end
   else
      error('Wrong mode!');
   end
end % function Aprod

function z = MIprod(x,mode)
    if mode == 1
        z1 = Aprod(W.*x,1);
        z  = [z1; damp*x];
    else
        z  = W.*Aprod(x(1:M),2) + damp*x(M+1:end);
    end
end

function z = Mprod(x,mode)
   if mode == 1
       z = Aprod(W.*x,1);
   elseif mode == 2
       z = W.*Aprod(x,2);
   else
      error('Wrong mode!');
   end
end % function Mprod

end % function SolveIRLS

