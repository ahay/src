%function y=medianfilter(x,n) 
%  
% The function MEDIANFILTER returns data smoothed by n-point median filter.
%
% Calling sequence-
% y=medianfilter(x,n)
%
% Input-
%	x	- column vector of input data x(n)
% 	n	- number of smoothing points
% Output-
%	y	- column vector of filtered data y(n)
%

function y=medianfilter(x,n) 

%----- Define the smoothing range indexes
lbound = floor(n/2);
rbound = n-lbound-1;

y=x;

%----- Run the filter
for i=lbound+1:length(x)-rbound-1
    y(i)=median(x(i-lbound:i+rbound));
end
