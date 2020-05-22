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
%
% Key Reference:
% Chen, et al., 2014,  Iterative deblending of simultaneous-source seismic data using seislet-domain shaping regularization, Geophysics, 79, V183-V193.
% Chen, Y., 2017, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 209, 21-31.
%
% Other related references (e.g., introducing the subroutines)
% Siahsar, M. A. N., Gholtashi, S., Kahoo, A. R., W. Chen, and Y. Chen, 2017, Data-driven multi-task sparse dictionary learning for noise attenuation of 3D seismic data, Geophysics, 82, V385-V396.
% Siahsar, M. A. N., V. Abolghasemi, and Y. Chen, 2017, Simultaneous denoising and interpolation of 2D seismic data using data-driven non-negative dictionary learning, Signal Processing, 141, 309-321.
% Chen, Y., M. Zhang, M. Bai, and W. Chen, 2019, Improving the signal-to-noise ratio of seismological datasets by unsupervised machine learning, Seismological Research Letters, 90, 1552-1564.
% Chen, Y., S. Zu, W. Chen, M. Zhang, and Z. Guan, 2019, Learning the blending spikes using sparse dictionaries, Geophysical Journal International, 218, 1379?1397. 
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51?KS61.
% Chen, Y., S. Fomel, 2015, Random noise attenuation using local signal-and-noise orthogonalization, Geophysics, 80, WD1-WD9.
% Chen, Y., J. Ma, and S. Fomel, 2016, Double-sparsity dictionary for seismic noise attenuation, Geophysics, 81, V17-V30.
% Zu, S., H. Zhou, R. Wu, M. Jiang, and Y. Chen, 2019, Dictionary learning based on dip patch selection training for random noise attenuation, Geophysics, 84, V169?V183.
% Zu, S., H. Zhou, R. Wu, and Y. Chen, 2019, Hybrid-sparsity constrained dictionary learning for iterative deblending of extremely noisy simultaneous-source data, IEEE Transactions on Geoscience and Remote Sensing, 57, 2249-2262.
% etc. 
 
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
