function [ dout ] = localsvd( d, param )
% Local SVD after dip steering 
% 
% d:            input data
% dout:         output data
% param:        parameters set
% param.r:      number of ranks in the shifted domain
% param.mode:   method in selecting the reference trace
%               1:      left trace
%               other:  averaged trace
%
% Reference: 
% Bekara, M., and M. van der Baan, 2007, Local singular value decomposition for signal
%   enhancement of seismic data: Geophysics, 72, V59–V65.
% 
% 
r=param.r; 
mode=param.mode;

if mode==1
    ref=d(:,1);
else
    nx=size(d,2);
    ref=sum(d,2)/nx;
end

  shift=dip_steering(d,ref);
  %figure;
  %plot(shift);
  d1=seisdither(d,shift);
  %figure;imagesc([d,d1]);
  
  %% inverse dip steering
  %d2=dither(d1,-shift);
  %figure;imagesc([d,d1,d2]);
  
  %% SVD after dip steering
  [u,e,v]=svd(d1);
  r=1;
  d1_svd=u(:,1:r)*e(1:r,1:r)*v(:,1:r)';
  dout=seisdither(d1_svd,-shift);

end

