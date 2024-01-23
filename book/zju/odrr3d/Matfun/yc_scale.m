function [ D1 ] = yc_scale(D,N,dscale)
%yc_scale: Scale the data up to the Nth dimension = sfscale axis=N
%  IN   D:   	intput data
%       N:      number of dimension for scaling
%               default: N=2
% 		dscale:  Scale by this factor
%       (does not include the rscale and pclip functions (not convenient actually))
%
%  OUT   D1:  	output data
%
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
%  Modified by Yangkang Chen on Jan, 2020
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%

if nargin==0
    error('Input data must be provided!');
end

if nargin==1
    N=2;
    dscale=1.0;
end

if nargin==2
    dscale=1.0;
end


[n1,n2,n3]=size(D);
D1=D;
switch N
    case 1
        for i3=1:n3
            for i2=1:n2
                D1(:,i2,i3)=D1(:,i2,i3)/max(abs(D1(:,i2,i3)));
            end
        end
    case 2
        for i3=1:n3
            D1(:,:,i3)=D1(:,:,i3)/max(max(abs(D1(:,:,i3))));
        end
    case 3
        D1=D1/max(max(max(abs(D1))));
        
    case 0
        
        D1=D1*dscale;
        
    otherwise
        error('Invalid argument value N.')
end


return
