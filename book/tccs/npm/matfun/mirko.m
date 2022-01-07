function mirko(Data_In, Data_Out,dt)

%=========================================================================
% Author      : Guoning Wu
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Jan, 2016

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
        
%  Copyright (C) 2016 The University of Texas at Austin
%  Copyright (C) 2016 Guoning Wu
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
%===========================================================================

dims = rsf_dim(Data_In);%demension of the input data
dd = zeros(dims(1), dims(2));
rsf_read(dd, Data_In);
imfs = dd';
[tf, tt, ff] = nnspe(imfs, -1.0, 10.23, 401,1124, 0, 40, -1, 10.23);

[tt1, ff1] = meshgrid(tt,ff);
tf1 = interp2(tt, ff, tf, tt1, ff1,'cubic');

q = fspecial('gaussian', 15, 5.0);
tf2 = filter2(q, tf1);

rsf_create(Data_Out, size(tf2)');
rsf_write(tf2.^.3, Data_Out);
