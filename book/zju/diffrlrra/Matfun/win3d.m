function [ dout ] =win3d(oper, param, din, n1win, n2win, n3win, r1, r2, r3)
% Processing in 3D windows
% 2D version similar to process_win
% coding strategy follows exactly usr/cm/Mfxmssa_win.c
%
% din:          input data
% oper:         operator
% param:        parameters of operator
% n1win:        first window length
% n2win:        second window length
% n3win:        third window length
% r1:           first overlapping ratio
% r2:           second overlapping ratio
% r3:           third overlapping ratio
%
% dout:     output data
% 
% Author      : Yangkang Chen
%
% Date        : Jan, 2018
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
% see also
% win2d.m,win3d_mask.m
%
% Example:
% test/test_win2d_fxmssa.m
% test/test_win3d_fxymssa.m
%
%  REFERENCES
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.



[n1,n2,n3]=size(din);

if nargin==3
    n1win=n1;
    n2win=n2;
    n3win=n3;
    r1=0.5;
    r2=0.5;
    r3=0.5;
end

if nargin==5
    r1=0.5;
    r2=0.5;
end

nov1=(1-r1)*n1win;  % non-overlapping size 1
nov2=(1-r2)*n2win;  % non-overlapping size 2
nov3=(1-r3)*n3win;  % non-overlapping size 3
ov1=r1*n1win;       % overlapping size 1
ov2=r2*n2win;       % overlapping size 2
ov3=r3*n3win;       % overlapping size 3

n1pad=n1win;        %padding size 1
nw1=1;
while n1pad<n1
    n1pad=n1pad+nov1;nw1=nw1+1;
end

n2pad=n2win;        %padding size 2
nw2=1;
while n2pad<n2
    n2pad=n2pad+nov2;nw2=nw2+1;
end

n3pad=n3win;        %padding size 3
nw3=1;
while n3pad<n3
    n3pad=n3pad+nov3;nw3=nw3+1;
end

D1=zeros(n1pad,n2pad, n3pad);D1(1:n1,1:n2,1:n3)=din;%copy din into padded D1
D2=zeros(n1pad,n2pad,n3pad);

for iw3=0:nw3-1
    for iw2=0:nw2-1
        for iw1=0:nw1-1
            s1=iw1*nov1;s2=iw2*nov2;s3=iw3*nov3;
            dtmp=D1(s1+1:s1+n1win,s2+1:s2+n2win,s3+1:s3+n3win);
            
            %uncomment this line for checking the correctness (checked 100% correct)
            dtmp = feval(oper,dtmp,param);
            %only valid for space-independent param
            %for reconstruction, with mask, param needs to be changed
            
            dtmp=win_weight3d(dtmp,iw1,iw2,iw3,nw1,nw2,nw3,n1win,n2win,n3win,ov1,ov2,ov3);
            
            D2(s1+1:s1+n1win,s2+1:s2+n2win,s3+1:s3+n3win)=D2(s1+1:s1+n1win,s2+1:s2+n2win,s3+1:s3+n3win)+dtmp;
        end
    end
end

dout=D2(1:n1,1:n2,1:n3);

return


function [dout ]= win_weight3d(din,iw1,iw2,iw3,nw1,nw2,nw3,n1win,n2win,n3win,ov1,ov2,ov3)
%       follow exactly usr/cm/win.c
%       float din /*input data*/,
% 		int iw1 /*starting window 1 in dst*/,
% 		int iw2 /*starting window 2 in dst*/,
% 		int iw3 /*starting window 3 in dst*/,
% 		int nw1 /*no of windows 1 in src*/,
% 		int nw2 /*no of windows 2 in src*/,
% 		int nw3 /*no of windows 3 in src*/,
% 		int n1win /*window length 1 in src*/,
% 		int n2win /*window legnth 2 in src*/,
% 		int n3win /*window legnth 3 in src*/,
% 		int ov1 /*copy length in axis1*/,
% 		int ov2 /*copy length in axis2*/,
% 		int ov3 /*copy length in axis3*/)


if iw3~=0
    for i1=0:n1win-1
        for i2=0:n2win-1
            for i3=0:ov3-1
                din(i1+1,i2+1,i3+1)=din(i1+1,i2+1,i3+1)*(i3+1)/(ov3+1);
            end
        end
    end
    
end

if iw3~=nw3-1
    for i1=0:n1win-1
        for i2=0:n2win-1
            for i3=0:ov3-1
                din(i1+1,i2+1,n3win-ov3+i3+1) = din(i1+1,i2+1,n3win-ov3+i3+1)*(ov3-i3)/(ov3+1);
            end
        end
    end
end



if iw2~=0
    for i3=0:n3win-1
        for i1=0:n1win-1
            for i2=0:ov2-1
                din(i1+1,i2+1,i3+1)=din(i1+1,i2+1,i3+1)*(i2+1)/(ov2+1);
            end
        end
    end
    
end

if iw2~=nw2-1
    for i3=0:n3win-1
        for i1=0:n1win-1
            for i2=0:ov2-1
                din(i1+1,n2win-ov2+i2+1,i3+1) = din(i1+1,n2win-ov2+i2+1,i3+1)*(ov2-i2)/(ov2+1);
            end
        end
    end
end

if iw1~=0
    for i3=0:n3win-1
        for i2=0:n2win-1
            for i1=0:ov1-1
                din(i1+1,i2+1,i3+1)=din(i1+1,i2+1,i3+1)*(i1+1)/(ov1+1);
            end
        end
    end
    
end

if iw1~=nw1-1
    for i3=0:n3win-1
        for i2=0:n2win-1
            for i1=0:ov1-1
                din(n1win-ov1+i1+1,i2+1,i3+1)=din(n1win-ov1+i1+1,i2+1,i3+1)*(ov1-i1)/(ov1+1);
            end
        end
    end
end

dout=din;
return












