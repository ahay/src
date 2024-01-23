function  [o] =  yc_bp(d,dt,f1,f2,f3,f4,verb);
%YC_BP: Apply a band-pass filter to a group of traces.
%
%  [o] = bp_filter(d,dt,f1,f2,f3,f4);
%
%  IN   d:    data (columns are traces)
%       dt:   sampling interval in sec
%       f1:   freq. in Hz
%       f2:   freq. in Hz
%       f3:   freq. in Hz
%       f4:   freq. in Hz
%
%   ^
%   |     ___________
%   |    /           \   Amplitude spectrum
%   |   /             \
%   |  /               \
%   |------------------------>
%      f1 f2        f3 f4
%
%  OUT  o:    output  (columns are traces)
%
%  Example: 
%
%    d=levents; 
%    dout = yc_bp(d,0.004,1,3,30,40); 
%    wigb([d,dout]);
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
% . 
%  The subroutine is simply renamed by Yangkang Chen, in Oct, 2019
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
[n1,n2,n3]=size(d);

if nargin==6
   verb=0; 
end

 [nt,nx] = size(d);
 k = nextpow2(nt);
 nf = 4*(2^k);

 i1 = floor(nf*f1*dt)+1;
 i2 = floor(nf*f2*dt)+1;
 i3 = floor(nf*f3*dt)+1;
 i4 = floor(nf*f4*dt)+1;

 up =  (1:1:(i2-i1))/(i2-i1);
 down = (i4-i3:-1:1)/(i4-i3);
 aux = [zeros(1,i1), up, ones(1,i3-i2), down, zeros(1,nf/2+1-i4) ];
 aux2 = fliplr(aux(1,2:nf/2));

 c = 0; % zero phase (could apply rotations as well)
 F = ([aux,aux2]');
 Phase = (pi/180.)*[0.,-c*ones(1,nf/2-1),0.,c*ones(1,nf/2-1)];
 
 Transfer = F.*exp(-i*Phase');


 D = fft(d,nf,1);

 for k = 1:nx
  Do(:,k) = Transfer.*D(:,k);
  
  if verb==1
     fprintf('ix/nx=%d/%d is done\n',k,nx); 
  end
 end

 o = ifft(Do,nf,1);

 o = real(o(1:nt,:));
 
 o=reshape(o,n1,n2,n3);


 
