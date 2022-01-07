 function Field(noisy,mssa,dmssa,n1,n2,n3,K,N)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Apr, 2015

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
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


%% from Madagascar to Matlab
data=zeros(n1,n2*n3);
rsf_read(data,noisy);

%% 2D to 3D
data=reshape(data,n1,n2,n3);

%% Noise attenuation
d1=fxymssa(data,0,120,0.004,K,0);	%MSSA
d2=fxydmssa(data,0,120,0.004,K,N,0);	%DMSSA


%% 3D to 2D
d1=reshape(d1,n1,n2*n3);
d2=reshape(d2,n1,n2*n3);


%% from Matlab to Madagascar
rsf_create(mssa,size(d1)');
rsf_write(d1,mssa);

rsf_create(dmssa,size(d2)');
rsf_write(d2,dmssa);



