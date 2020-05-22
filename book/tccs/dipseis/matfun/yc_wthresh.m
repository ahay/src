function y = yc_wthresh(x,sorh,t) 
%WTHRESH Perform soft or hard thresholding.  
%   Y = WTHRESH(X,SORH,T) returns soft (if SORH = 's') 
%   or hard (if SORH = 'h') T-thresholding  of the input  
%   vector or matrix X. T is the threshold value. 
% 
%   Y = WTHRESH(X,'s',T) returns Y = SIGN(X).(|X|-T)+, soft  
%   thresholding is shrinkage. 
% 
%   Y = WTHRESH(X,'h',T) returns Y = X.1_(|X|>T), hard 
%   thresholding is cruder. 
% 
%   See also WDEN, WDENCMP, WPDENCMP. 
 
%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96. 
%   Last Revision: 14-May-2003. 
%   Copyright 1995-2004 The MathWorks, Inc. 
% $Revision: 1.11.4.2 $ 
 
switch sorh 
  case 's' 
    tmp = (abs(x)-t); 
    tmp = (tmp+abs(tmp))/2; 
    y   = sign(x).*tmp; 
  
  case 'h' 
    y   = x.*(abs(x)>t); 
  
  otherwise 
    error('Invalid argument value.') 
end
