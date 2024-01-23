function FAULTS(din,dout)
% Author      : Hang Wang, Xingye Liu, and Yangkang Chen
% Date        : March, 2020

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API

%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Hang Wang, Xingye Liu, and Yangkang Chen
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
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
%
%  REFERENCES
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.


d=zeros(1000,250);
rsf_read(d,din);

d1=fxymssa_win_auto2(d,0,120,0.004,3,6,0,200,50,1,0.5,0.5,0.5,2);

rsf_create(dout,size(d1)');
rsf_write(d1,dout);


    