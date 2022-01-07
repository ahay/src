function [ D1 ] = yc_fkt(D,sorh,t)
%YC_FKT: FK domain thrsholding (any dimension)
%  IN   D:   	intput data 
%       sorh:   soft or hard thresholding (s,h,ps,ph)
%               s: soft
%               h: hard
%               ps: percentile soft
%               ph: percentile hard
%       t:      threshold value
%      
%  OUT   D1:  	output data
% 
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
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
if nargin==1 % for quick denoising test
   sorh='ps';
   t=5;
end
D1=ifftn(yc_pthresh(fftn(D),sorh,t));

% [n1,n2,n3]=size(D);
% n1h=floor(n1/2)+1;
% 
% Dfft=reshape(fftn(D),n1,n2*n3);
% Dtmp=reshape(Dfft(1:n1h,:),n1h,n2,n3); % thresholding only on the half space
% %Dfft(1:n1h,:)=reshape(pthresh(Dtmp,sorh,t),n1h,n2*n3);
% 
% Dfft(1:n1h,:)=reshape(Dtmp,n1h,n2*n3);
% 
% %% honor conjugate symmetry property
% if mod(n1,2)==0 % even
% Dfft(n1h+1:end,:)=conj(flipud(Dfft(2:n1h-1,:)));
% else            % odd
% Dfft(n1h+1:end,:)=conj(flipud(Dfft(2:n1h,:)));    
% end
% Dfft=reshape(Dfft,n1,n2,n3);
% %D1=ifftn(Dfft);
% D1=Dfft;
return
