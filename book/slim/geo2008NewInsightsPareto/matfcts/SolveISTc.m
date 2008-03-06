function [x, xnorms, rnorms, lambdas] = SolveISTc(A, y, Iters, InnerIters, x, fullPath)
% SolveISTc: Iterative Soft Thresholding with cooling
%
%---------------------------------------------------------------
% Solve the basis pursuit (BP) problem
%
%      min ||x||_1  s.t.  Ax=y
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
% b             is an m-vector.
% Iters         Total number of iterations
% InnerIters    Number of iterations per subproblem
% x             is an N-vector estimate of the solution (possibly all
%               zeros). If x = [], then ISTc determines N via
%               N = length( A'b ) and sets  x = zeros(N,1).
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

%----------------------------------------------------------------------
% Check arguments. 
%----------------------------------------------------------------------
if nargin < 6
    fullPath = 0;
end

if nargin < 5
    x = [];
end

if nargin < 4
    InnerIters = 5;
end

if nargin < 3
    Iters = 20;
end

if nargin < 2
    error('At least two arguments are required');
end

OuterIters = ceil(Iters/InnerIters);
iter = 0;

xnorms  = [];
rnorms  = [];
lambdas = [];

%----------------------------------------------------------------------
% Determine initial x and vector length N
%----------------------------------------------------------------------
explicit = isnumeric(A);
object   = isobject(A);

if explicit || object
    projres = A'*y;
else
    projres = Aprod(y,2);
end
N = length(projres);

if isempty(x)
    x = zeros(N,1);
end

rnorms(1)  = norm(y-Aprod(x,1),2);
xnorms(1)  = norm(x,1);
lambdas(1) = 0;

% Compute the threshold values for the subproblems
maxthr = .99*max(abs(projres));

%----------------------------------------------------------------------
% MAIN LOOP
%----------------------------------------------------------------------
h = waitbar(0,'SolveISTc - Please wait...');
for i = 1:OuterIters
    thr = maxthr/1.5^(i-1); % popular cooling strategy
    for j = 1:InnerIters
        xTmp = SoftThr(x + projres - Aprod(Aprod(x,1),2),thr);
        if fullPath
            rnorms((i-1)*InnerIters+j+1) = norm(y-Aprod(xTmp,1),2);
            xnorms((i-1)*InnerIters+j+1) = norm(xTmp,1);
            lambdas((i-1)*InnerIters+j+1) = thr;
        end
        x = xTmp;
        iter = iter+1;
        waitbar(iter/Iters,h)
        if iter == Iters
            close(h)
            break
        end
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

end % function SolveISTc

function  x = SoftThr(y,thr)
    ay = abs(y);
    x = sign(y).*(ay >= thr).*(ay - thr);
end % SoftThr
