function Curve(infile,N1,N2,outfile1,outfile2)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%
% Date        : Feb, 2015

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

%% infile: 	input seismic image
%% outfile1: 	coeficients
%% outfile2: 	display of coefficients

%% default values
    X=zeros(N1,N2);
    rsf_read(X,infile);
    
    C = yc_fdct_wrapping(X,0,2); 

   %generate curvelet image (a complex array)
    img = yc_fdct_wrapping_dispcoef(C);
    img=abs(img);

    nscale=max(size(C));
    
    coef=[];
    for iscale=1:nscale
    	   nangle=max(size(C(iscale)));
		for iangle=1:nangle
			[n1,n2]=size(C{iscale}{iangle});
			temp=reshape(C{iscale}{iangle},n1*n2,1);
			coef=[coef;temp];
    		end
    end		
    coef=abs(coef);
          
    rsf_create(outfile1,size(coef)');
    rsf_write(coef,outfile1);

    %% from Matlab to Madagascar
    rsf_create(outfile2,size(img)');
    rsf_write(img,outfile2);







