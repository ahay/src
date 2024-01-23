 function FXY_MSSA_WIN_AUTO(din,d_lrra,n1,n2,n3,dt,lf,hf,N,n1win,n2win,n3win,r1,r2,r3,verb)
% Author      : Yangkang Chen
%         
% Date        : Nov, 2020

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2020 Zhejiang University
%  Copyright (C) 2020 Yangkang Chen
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
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.

%% from Madagascar to Matlab
d=zeros(n1,n2*n3);
rsf_read(d,din)
d=reshape(d,n1,n2,n3);

%%% Main program goes here !
    %% LDRR
    d2=fxymssa_win_auto(d,lf,hf,dt,N,verb,n1win,n2win,n3win,r1,r2,r3,2);
    d2=reshape(d2,n1,n2,n3);

%% from Matlab to Madagascar
    rsf_create(d_lrra,size(d2)');
    rsf_write(d2,d_lrra);

return

