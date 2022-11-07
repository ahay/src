

function VectaTraceNar(DataInput, TimeFreqCube, TimeFreqSlice1, TimeFreqSlice2,...
                    TimeFreqSlice3, TimeFreqSlice4, TimeFreqSlice5, TimeFreqSlice6)

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

% parameters setting

nt = 751;
dt = 0.002;
t0 = 0.0;
t1 = (nt-1)*dt;
timeres = nt;

fs = 1/dt;
nf = 128;
df = (fs/2/nf);
f1 = (nf-1)*df;
f0 = 0;
freres = fs;

% from madagascar --> matlab

dims = rsf_dim(DataInput);             %demension of the input data
mat_data = zeros(dims(1), dims(2));
rsf_read(mat_data, DataInput);

% matlab process

slice20 = [];
slice40 = [];
slice60 = [];
slice80 = [];
slice100 = [];
slice120 = [];
tfemd = [];
tic;
for i=1:5:dims(2)-2
    tf = [];
    imfs = mat_data(:,i:i+4);
    [tf, tt, ff] = nnspe(imfs, t0, t1, freres, timeres, f0, f1, t0, t1);
    [tt1,ff1] = meshgrid(tt,ff);

    tf = interp2(tt,ff,tf,tt1,ff1,'cubic');
    q = fspecial('gaussian', 15, 3.0);
    tf = filter2(q, tf);
    tf = filter2(q, tf);
    tfemd = [tfemd,tf'];
    slice20 = [slice20, tf(20,:)'];
    slice40 = [slice40, tf(40,:)'];
    slice60 = [slice60, tf(60,:)'];
    slice80 = [slice80, tf(80,:)'];
    slice100 = [slice100, tf(100,:)'];
    slice120 = [slice120, tf(120,:)'];
end
toc;
% from MatLab --> Madagascar

rsf_create(TimeFreqCube, size(tfemd)');
rsf_write(tfemd, TimeFreqCube);

rsf_create(TimeFreqSlice1, size(slice20)');
rsf_write(slice20, TimeFreqSlice1);

rsf_create(TimeFreqSlice2, size(slice40)');
rsf_write(slice40, TimeFreqSlice2);

rsf_create(TimeFreqSlice3, size(slice60)');
rsf_write(slice60, TimeFreqSlice3);

rsf_create(TimeFreqSlice4, size(slice80)');
rsf_write(slice80, TimeFreqSlice4);

rsf_create(TimeFreqSlice5, size(slice100)');
rsf_write(slice100, TimeFreqSlice5);

rsf_create(TimeFreqSlice6, size(slice120)');
rsf_write(slice120, TimeFreqSlice6);

end
