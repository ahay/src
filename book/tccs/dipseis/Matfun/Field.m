function Field(data,planes)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Apr, 2014

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
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


% clear;clc;close all;

dt=0.001;
verb=1;
d1=zeros(700,310);
rsf_read(d1,data)



dip1=d1-fxemd(d1,5, 450, 0.001, 1, verb);
dip2=d1-fxemd(d1,5, 450, 0.001, 2, verb)-dip1;
dip3=d1-fxemd(d1,5, 450, 0.001, 3, verb)-dip1-dip2;
dip4=d1-fxemd(d1,5, 450, 0.001, 4, verb)-dip1-dip2-dip3;
res=d1-dip1-dip2-dip3-dip4;
%figure;imagesc([[d1,dip1,dip2];[dip3,dip4,res]]);

dips=[d1,dip1,dip2,dip3,dip4,res];

%rsf_create('dips.rsf',size(dips)')
%rsf_write(dips,'dips.rsf');

rsf_create(planes,size(dips)');
rsf_write(dips,planes);


%figure;imagesc([[d1,dip1,dip2];[dip3,dip4,res]]);
% plane1=dip1+dip2;
% plane2=dip3+dip4+res;
% figure;imagesc([plane1,plane2]);

%figure;imagesc([plane1,plane2,plane3]);

%figure;imagesc([d1,dip1,dip2;[dip3,dip4,res]]);
