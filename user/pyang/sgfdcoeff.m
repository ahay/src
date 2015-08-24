function c=sgfdcoeff(NJ)
% compute 2*N-order staggered-grid finite difference coefficients (NJ=2N)
% Example:
%     format long
%     NJ=10;
%     c=staggered_fdcoeff(NJ)
% 
% Note: This code is devoted alongside the manuscript of P.L. Yang et al, 
% Using the effective boundary saving strategy in GPU-based RTM programming
%
% Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)

N=NJ/2;
x=zeros(N,1);
b=zeros(N,1);   b(1)=1;
c=b;
for k=1:N
    x(k)=(2*k-1)^2;
end

for k=1:N-1
    for i=N:-1:k+1
        b(i)=b(i)-x(k)*b(i-1);
    end
end

for k=N-1:-1:1
    for i=k+1:N
        b(i)=b(i)/(x(i)-x(i-k));
    end
    for i=k:N-1
        b(i)=b(i)-b(i+1);
    end
end
for k=1:N
    c(k)=b(k)/(2*k-1);
end
