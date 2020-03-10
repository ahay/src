function [ shift ] = dip_steering( din, ref)
% Dip steering
%   Detailed explanation goes here

%  din: input data
%  ref: reference trace
%  dout: output data

[nt,nx]=size(din);
nlag=300;
shift=zeros(1,nx);
    for ix=1:nx
           [corr,lag]=xcorr(ref,din(:,ix),nlag);
           [vmax,lagmax]=max(abs(corr));
           shift(ix)=lag(lagmax);        
    end
end

