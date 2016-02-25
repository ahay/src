function [dout] = ddtf_denoise2d(din, lambda, niter, lvl, htype, thr, thrtype, sh)
%DDTF_DENOISE2D: Data-driven tight frame denoising for 2D seismic data
%
%	 [ D1 ] = ddtf_denoise2d(din, lambda, niter, lvl, htype, thr, thrtype);
%
%  IN       din:        input data
%           lambda:     controling parameter
%           niter:      linearing iteration
%           lvl:        filter size (2^lvl-1)
%           htype:      filter type 
%           thr:        thresholding level
%			thrtype:	threshold type (g-global, l-local, e-exact)
% 	    sh: 	soft or hard
%               
%  OUT      dout:  	output data
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
if nargin==0
 error('Input data must be provided!');
end

if nargin==1
    lambda=0.5;
    niter=30;
    lvl=2;          %patch size=7*7
    htype='spline';
    thr=20;       %input data should better be scaled
	thrtype='g';	
end;

if nargin==7
   sh='s';
end

% initial basis
[H0, hsize] = initialFilter(htype, lvl);



% learning filters
H = ddTF(din, H0,hsize, lambda,niter);
d = ddtf_dec_p(din, H);

%% 
if thrtype=='e'
    dthr = ani_thresh(d, sh, thr);
else
    dthr = ani_thresh_cyk(d, sh, thr, thrtype);
end

dout=ddtf_rec_p(dthr,H);
   
end
