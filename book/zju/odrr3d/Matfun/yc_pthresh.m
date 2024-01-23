function [y,thr] = yc_pthresh(x,sorh,t) 
%YC_PTHRESH Perform soft or hard thresholding or percentile
%        soft or hard thresholding.  
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
%
%   Yangkang Chen, The University of Texas at Austin

 
switch sorh 
  case 's' 
    tmp = (abs(x)-t); 
    tmp = (tmp+abs(tmp))/2; 
    y   = sign(x).*tmp; 
  
  case 'h' 
    y   = x.*(abs(x)>t);
    
  case 'ps'
   tmp=reshape(abs(x),1,prod(size(x)));
   t=prctile(tmp,100-t);       
   if nargout==2
      thr=t; 
   end
    tmp = (abs(x)-t); 
    tmp = (tmp+abs(tmp))/2; 
    y   = sign(x).*tmp; 
  case 'ph'    
   tmp=reshape(abs(x),1,prod(size(x)));
   t=prctile(tmp,100-t);       
   if nargout==2
      thr=t; 
   end
    y   = x.*(abs(x)>t);   
  otherwise 
    error('Invalid argument value.') 
end
