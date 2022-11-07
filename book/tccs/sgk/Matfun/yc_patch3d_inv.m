function [ A ] = yc_patch3d_inv( X,mode,n1,n2,n3,l1,l2,l3,s1,s2,s3 )
% insert patches into the 3D data
%
% by Yangkang Chen
% March, 2020
%
% Input:
%   D: input image
%   mode: patching mode
%   l1: first patch size
%   l2: second patch size
%   l3: third patch size
%   s1: first shifting size
%   s2: second shifting size
%   s3: third shifting size
%
% Output:
%   X: patches
%
% Modified on Dec 12, 2018 (the edge issue, arbitrary size for the matrix)
% 		      Dec 31, 2018 (tmp1=mod(n1,l1) -> tmp1=mod(n1-l1,s1))
%             Marich, 31, 2020, 2D->3D
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

%% patch size l1*l2*l3
%l1=4;l2=4;l3=4;
%

if mode==1 %possible for other patching options
    
    if nargin==4
        l1=4;l2=4;l3=4;s1=4;s2=4;s3=4;
    end
    
    if nargin==6
        s1=round(l1/2);s2=round(l2/2);s3=round(l3/2);
    end
    
    tmp1=mod(n1-l1,s1);
    tmp2=mod(n2-l2,s2);
    tmp3=mod(n3-l3,s3);
    if tmp1~=0 && tmp2~=0 && tmp3~=0
        A=zeros(n1+s1-tmp1,n2+s2-tmp2,n3+s3-tmp3);
        mask=zeros(n1+s1-tmp1,n2+s2-tmp2,n3+s3-tmp3);
    end
    
    if tmp1~=0 && tmp2==0 && tmp3==0
        A=zeros(n1+s1-tmp1,n2,n3);
        mask=zeros(n1+s1-tmp1,n2,n3);
    end
    
    if tmp1==0 && tmp2~=0 && tmp3==0
        A=zeros(n1,n2+s2-tmp2,n3);
        mask=zeros(n1,n2+s2-tmp2,n3);
    end
    
    if tmp1==0 && tmp2==0 && tmp3~=0
        A=zeros(n1,n2,n3+s3-tmp3);
        mask=zeros(n1,n2,n3+s3-tmp3);
    end
    
    if tmp1==0 && tmp2==0  && tmp2==0
        A=zeros(n1,n2,n3);
        mask=zeros(n1,n2,n3);
    end
    
    [N1,N2,N3]=size(A);
    id=0;
    for i1=1:s1:N1-l1+1
        for i2=1:s2:N2-l2+1
            for i3=1:s3:N3-l3+1
                id=id+1;
                A(i1:i1+l1-1,i2:i2+l2-1,i3:i3+l3-1)=A(i1:i1+l1-1,i2:i2+l2-1,i3:i3+l3-1)+reshape(X(:,id),l1,l2,l3);
                mask(i1:i1+l1-1,i2:i2+l2-1,i3:i3+l3-1)=mask(i1:i1+l1-1,i2:i2+l2-1,i3:i3+l3-1)+ones(l1,l2,l3);
            end
        end
    end
    
    A=A./mask;
    
    A=A(1:n1,1:n2,1:n3);
end

return

