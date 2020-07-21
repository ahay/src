%function [f, a] = FAimphilbert(data,dt)

% The function FAimphilbert calculates frequency and amplitude
% of data(n,k) using an improved Hilbert method, where n specifies the length
%  of time series, and k is the number of IMF components.
%    hilbt.m is used to perform Hilbert transform instead of hilbert.m,
%     which mainly reduce the Gipps impacts from the end points.
% 
% hilbtm.m is used to perform an improved Hilbert transform.(Gibbs phenomena is fixed)
% Note: FAH allows the instantaneous frequency to be negative. 
%
% Calling sequence-
% [f,a] = FAimphilbert(data,dt)
%
% Input-
%   data	- 2-D matrix data(n,k) of IMF components 
%	  dt	  - sampling period in seconds
% Output-
%	  f	    - 2-D matrix f(n,k) that specifies the Hilbert frequency in Hz
%	  a	    - 2-D matrix a(n,k) that specifies the Hilbert amplitude
%
% Used by-
% 	FA, NSPABMUN.
%
%written by  
% Norden Huang (NASA GSFC) Junw 2,2002 :Initial
% Xianyao Chen, september, 2008 : modified
% S.C.Su ,Sep ,2009, rename all the functions
%
function [f, a] = FAimphilbert(data,dt)

%----- Get the dimension
%[nPoints, nIMF] = size(data);
[nIMF,npt] = size(data); 
flip=0;  
if nIMF>npt
    data=data';
    [nIMF,npt] = size(data);
    flip=1;
end   

%----- Apply improved Hilbert transform
h=hilbtm(data);

%----- Get the instantaneous frequency
f(1:nIMF,1:npt)=0;
for j1=1:nIMF
temp=diff(unwrap(angle(h(j1,:))))./(2*pi*dt);

%----- Duplicate last points to make dimensions equal
tta=[temp temp(end)];
f(j1,1:npt)=tta(1:npt);
clear tta temp 
end

%----- Get the amplitude
a=abs(h);
if flip ==1
    f=f';
    a=a';
end    
