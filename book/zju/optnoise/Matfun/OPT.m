 function OPT(infile,dip,outfile,n1,n2,N,r)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Jun, 2016

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2016 The University of Texas at Austin
%  Copyright (C) 2016 Yangkang Chen
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

%%from Madagascar to Matlab
% create memory
d1=zeros(n1,n2);
rsf_read(d1,infile)

sigma=zeros(n1,n2);
rsf_read(sigma,dip)

%% main program
p=orthogonal_polynomial([1:n2],N);

% for adjusting the flattening error
d2=zeros(n1,n2);
for i2=1:n2
    d2(:,i2)=yc_pwtp(d1(:,i2),sigma,i2,1);  
end

% masking
%r=2;
cn=d2*p'*inv(p*p');
cnp=zeros(n1,N);
cnp(:,1:r)=cn(:,1:r);
d3=cnp*p;

%% reverse adjustment
dout=zeros(n1,n2);
for i2=1:n2
    dout(:,i2)=yc_pwtp(d3(:,i2),-sigma,i2,1);  
end

rsf_create(outfile,size(dout)');
rsf_write(dout,outfile);
