function gch1(data,tfsswt)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Nov, 2013

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

dt=0.008;
nt=501;
of=0;
nf=257;
fs=1/dt;
df=(fs/2-of)/(nf-1);

% allocate memory
d = zeros(1,nt);

% from Madagascar to Matlab
rsf_read(d,data);

[dtemp, f0, opts]=sswt(d, fs,'Display','notify');
%[d2, f, opts]=sswft(d, 125,'Display','notify');
%[d3, f, opts]=wt(d, 125,'Display','notify');
%[d4, f, opts]=wft(d, 125,'Display','notify');
t0=[0:500]*dt;
f1=[0:nf-1]*df;

[t1,f1]=meshgrid(t0,f1);
dtf=interp2(t0,f0,dtemp,t1,f1);

rsf_create(tfsswt,size(dtf)');
rsf_write(abs(dtf),tfsswt);

















