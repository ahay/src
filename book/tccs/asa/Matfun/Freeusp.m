function Freeusp(Din,n1,n2,D1,D2,D3,D4)
% Author      : Yangkang Chen
%				Zhejiang University
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
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

din=zeros(n1,n2);
d1=zeros(n1,n2);
d2=zeros(n1,n2);
d3=zeros(n1,n2);
d4=zeros(n1,n2);
rsf_read(din,Din);

% denoise
for i2=1:n2
%	temp=emd(din(:,i2));
	temp=emd(din(:,i2),'MAXMODES',5);
	d1(:,i2)=temp(1,:);
	d2(:,i2)=temp(2,:);
	d3(:,i2)=temp(3,:);
	d4(:,i2)=temp(4,:);
	fprintf('Trace %d /%d is done !\n\n',i2,n2);
end

% from Matlab to Madagascar
rsf_create(D1,size(d1)');
rsf_write(d1,D1);

rsf_create(D2,size(d2)');
rsf_write(d2,D2);

rsf_create(D3,size(d3)');
rsf_write(d3,D3);

rsf_create(D4,size(d4)');
rsf_write(d4,D4);



















