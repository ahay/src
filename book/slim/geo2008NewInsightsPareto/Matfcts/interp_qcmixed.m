function yi = interp_qcmixed(x,f,g,xi,alpha)
% Authors     : E. van den Berg and M. P. Friedlander
%               Scientific Computing Laboratory
%               Department of Computer Science
%               The University of British Columbia
%         
% Date        : March, 08

%  Copyright (C) 2007 The University of British Columbia at Vancouver
%  Copyright (C) 2007 E. van den Berg and M. P. Friedlander
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

    yi = zeros(size(xi));

   % ---------------------------------------------------------

   % Quadratic interpolation - solve system for coefficients
   for i=1:length(x)-1
     x1 = x(i);
     x2 = x(i+1);
     A  = [1, x1, x1^2 ; ...
           1, x2, x2^2 ; ...
           0, 1,  2*x2];
     b  = [f(i); f(i+1); g(i+1)];
     c  = pinv(A)*b;
     idx = find((xi >= x1) & (xi <= x2));
     if ~isempty(idx)
        xs  = xi(idx);
        ys  = c(1) + c(2) * xs + c(3) * xs.^2;
        yi(idx) = ys;
     end
   end
      
   % ---------------------------------------------------------
   
   % Cubic interpolation - solve system for coefficients
   for i=1:length(x)-1
     x1 = x(i);
     x2 = x(i+1);
     A  = [1, x1, x1^2, x1^3  ; ...
           0, 1,  2*x1, 3*x1^2; ...
           1, x2, x2^2, x2^3  ; ...
           0, 1,  2*x2, 3*x2^2];
     b  = [f(i); g(i); f(i+1); g(i+1)];
     c  = pinv(A)*b;

     % Check convexity
     xr = -2*c(3) / (2 * 3 * c(4));
     if (xr >= 1.04*x1) && (xr <= 0.96 * x2)
        % Significant over-fit detected
        fprintf('Segment %d; restoring quadratic fit\n',i);
        continue; % Keep quadratic interpolation
     end
     

     idx = find((xi >= x1) & (xi <= x2));
     if ~isempty(idx)
        xs  = xi(idx);
        ys  = c(1) + c(2) * xs + c(3) * xs.^2 + c(4) * xs.^3;
        yi(idx) = ys;
     end
   end

   % ---------------------------------------------------------
   
   % Linear extrapolation
   al  = 0;
   cl  = f(end) - x(end) * g(end);
   bl  = g(end);
   
   % Quadratic extrapolation
   a = g(end).^2 / (4 * f(end));
   b = g(end) - 2 * x(end) * a;
   c = f(end) - x(end) * g(end) + x(end).^2 * a;
   
   % Linear combination of linear and quadratic
   a = (alpha * a + (1 - alpha) * al); 
   b = (alpha * b + (1 - alpha) * bl);
   c = (alpha * c + (1 - alpha) * cl);

   idx = find(xi >= x(end));
   xs  = xi(idx);
   ys = a * xs.^2 + b * xs + c;
   yi(idx) = ys;
      
end
