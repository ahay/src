function [ X ] = yc_patch( A,mode,l1,l2,s1,s2 )
%decompose the image into patches:
%  
% by Yangkang Chen
% Oct, 2017
%
% Input: 
%   D: input image
%   mode: patching mode
%   l1: first patch size
%   l2: second patch size
%   s1: first shifting size
%   s2: second shifting size
%   
% Output:
%   X: patches
% 
% Modified on Dec 12, 2018 (the edge issue, arbitrary size for the matrix)
% 		      Dec 31, 2018 (tmp1=mod(n1,l1) -> tmp1=mod(n1-l1,s1))
% 
%
% Key Reference:
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

%% patch size l1*l2
%l1=8;l2=8;
%

[n1,n2]=size(A);

if mode==1 %possible for other patching options

if nargin==2
   l1=8;l2=8;s1=4;s2=4; 
end

if nargin==4
   s1=round(l1/2);s2=round(l2/2);
end

tmp=mod(n1-l1,s1);
if tmp~=0
   A=[A;zeros(s1-tmp,n2)]; 
end
tmp=mod(n2-l2,s2);
if tmp~=0
   A=[A,zeros(size(A,1),s2-tmp)]; 
end



[N1,N2]=size(A);
 X=[];
for i1=1:s1:N1-l1+1
    for i2=1:s2:N2-l2+1
        tmp=reshape(A(i1:i1+l1-1,i2:i2+l2-1),l1*l2,1);
        X=[X,tmp];  
    end
end   
    
end




end

