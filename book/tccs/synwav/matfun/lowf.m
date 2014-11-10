function lowf(data,tf1,tf2,tf3,tf4,tf5,tf6,tf)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Dec, 2013

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2013 The University of Texas at Austin
%  Copyright (C) 2013 Yangkang Chen
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

%dt=0.004;
%nt=501;
%nf=126;

dt=0.002;
nt=751;
dx=1;
nx=201;

nf=126

fs=1/dt;
df=fs/2/(nf-1);

% allocate memory
d = zeros(nt,nx);

% from Madagascar to Matlab
rsf_read(d,data);

d1=[];
d2=[];
d3=[];
d4=[];
d5=[];
d6=[];
dd=[];

for i=1:nx
	[dtemp, f0, opts]=sswt(d(:,i), fs,'Display','notify');
	t0=[0:nt-1]*dt;
	f1=[0:nf-1]*df;
	[t1,f1]=meshgrid(t0,f1);
	dtf=interp2(t0,f0,dtemp,t1,f1);
	d1=[d1,abs(dtf(20,:))'];
	d2=[d2,abs(dtf(40,:))'];
	d3=[d3,abs(dtf(60,:))'];
	d4=[d4,abs(dtf(80,:))'];
	d5=[d5,abs(dtf(100,:))'];
	d6=[d6,abs(dtf(120,:))'];
	dd=[dd,abs(dtf)'];
	fprintf('Trace %d is done \n\n',i);
end

rsf_create(tf1,size(d1)');
rsf_write(d1,tf1);

rsf_create(tf2,size(d2)');
rsf_write(d2,tf2);

rsf_create(tf3,size(d3)');
rsf_write(d3,tf3);

rsf_create(tf4,size(d4)');
rsf_write(d4,tf4);

rsf_create(tf5,size(d5)');
rsf_write(d5,tf5);

rsf_create(tf6,size(d6)');
rsf_write(d6,tf6);

rsf_create(tf,size(dd)');
rsf_write(dd,tf);
















