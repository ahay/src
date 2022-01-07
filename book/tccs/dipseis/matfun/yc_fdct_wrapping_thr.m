 function [dout] = yc_fdct_wrapping_thr(din,perc);
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Mar, 2015

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


[n1,n2]=size(din);

C = yc_fdct_wrapping(din,1,1);
C1=C;
cv=[];%vector form of C
l1=length(C);
for i1=1:l1
    l2=length(C{i1});	
    for i2=1:l2
	cv=[cv;C{i1}{i2}(:)];
    end	
end
thr=prctile(cv,100-perc); 

for i1=1:l1	
    l2=length(C{i1});	
    for i2=1:l2
	C1{i1}{i2}=yc_wthresh(C{i1}{i2},'s',thr);
    end	
end

dout = real(yc_ifdct_wrapping(C1,1,n1,n2)); 
