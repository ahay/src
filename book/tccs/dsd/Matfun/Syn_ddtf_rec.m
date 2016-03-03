function Syn_ddtf_rec(din,din2,dout,n1,n2,lambda,niter,lvl,htype)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : May, 2014

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
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  


n23=n2*(2^(lvl+1)-1)^2;
un=zeros(n1,n23); % noisy input
rsf_read(un,din);

un2=zeros(n1,n2); % noisy input
rsf_read(un2,din2);

%if nargin==2
%    lambda=0.5;
%    niter=30;
%    lvl=2;          %patch size=7*7
%    htype='spline';
%    thr=20;       %input data should better be scaled
%	thrtype='g';	
%end;

% initial basis
[H0, hsize] = initialFilter(htype, lvl);

% learning filters
H = ddTF(un2, H0,hsize, lambda,niter);

% 
m=1;
n=n23/n2;
d = cell(m,n);

for i=1:n
	d{i}=un(:,i+(n2-1)*(i-1):i*n2);
end
u1=ddtf_rec_p(d,H);

% from Matlab to Madagascar
rsf_create(dout,size(u1)');
rsf_write(u1,dout);




















