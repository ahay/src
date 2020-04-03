 function Gendata(clean,noisy)
% Simulating AVO data
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Jun, 2015

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

%% make data (AVO)
eps=0.000000001;
dt=0.004;
nt=151;
dx=25;
nx=61;
tau1=0.1;p1=0;
tau2=0.2;p2=0;
tau3=0.3;p3=0;
tau4=0.4;p4=0;
tau5=0.5;p5=0;


avo=1-[1:nx]*2/(nx-1);
avo2=2-[1:nx]/(nx-1);
avo3=0.5+0.5*[1:nx]/(nx-1);
avo4=2-[1:nx]/(nx-1);
%plot(avo);
d=zeros(nt,nx);
for it=10:nt
    for ix=1:nx
        t=(it-1)*dt;
        x=(ix-1)*dx;
        if(abs(t-(tau1+p1*x^2))<eps)
           
           d(it,ix)=d(it,ix)+ avo(ix);
        end
        if(abs(t-(tau2+p2*x^2))<eps)
           
           d(it,ix)=d(it,ix)+ avo2(ix);
        end
        if(abs(t-(tau3+p3*x^2))<eps)
           
           d(it,ix)=d(it,ix)+ avo3(ix);
                end
        if(abs(t-(tau4+p4*x^2))<eps)
           
           d(it,ix)=d(it,ix)+ avo4(ix);
        end
        if(abs(t-(tau5+p5*x^2))<eps)
           
           d(it,ix)=d(it,ix)+ avo(ix);
        end

    end
end

w = yc_ricker(40,dt);
d = conv2(d,w,'same');
d=d/max(max(d));

randn('state',201314')
dn=d+0.2*randn(nt,nx);


%% from Matlab to Madagascar
rsf_create(clean,size(d)');
rsf_write(d,clean);

rsf_create(noisy,size(dn)');
rsf_write(dn,noisy);




