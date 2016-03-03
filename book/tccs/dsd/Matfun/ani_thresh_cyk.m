function [d,coef] = ani_thresh_cyk(c, sorh, thr, thrtype)
%ANI_THRESH_CYK: A modified version of ani_thresh
%
%	 [ D1 ] = ani_thresh_cyk(c, sorh, thr, thrtype)
%
%  IN       c:		input coefficients
%			sorh:   'soft' or 'hard'
%			thr:	percentile threshold value
%			thrtype:percentile threshold type
%              
%  OUT      d:  	output coefficients
%
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
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

[m, n] = size(c);
[m0,n0]=size(c{1});

d = cell(m,n);

coef=[];
for k=1:m*n
    coef=[coef,c{k}];
end
   
if thrtype=='g'
   temp=reshape(coef,1,m*n*m0*n0);
   lambda=prctile(temp,100-thr);     
   d{1} = c{1};
    for k = 2:m*n
        d{k} = wthresh(c{k}, sorh, lambda);
    end
else     
   d{1} = c{1};
    for k = 2:m*n
        temp=reshape(c{k},1,m0*n0);   
        lambda=prctile(temp,100-thr);  
        d{k} = wthresh(c{k}, sorh, lambda);
    end        
end

end

