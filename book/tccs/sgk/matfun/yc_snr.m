function psnr = yc_snr(g,f,mode)
% Author: Yangkang Chen
% g: ground truth image
% f: noisy/restored image
% mode:1->2D SNR, 2->3D SNR
% 
% References:
% Chen, Y., S. Fomel, 2015, Random noise attenuation using local signal-and-noise orthogonalization, Geophysics, 80, WD1-WD9.
% Chen, Y., J. Ma, and S. Fomel, 2016, Double-sparsity dictionary for seismic noise attenuation, Geophysics, 81, V17-V30.
% Chen, Y., 2017, Fast dictionary learning for noise attenuation of multidimensional seismic data, Geophysical Journal International, 209, 21-31.
% Siahsar, M. A. N., Gholtashi, S., Kahoo, A. R., W. Chen, and Y. Chen, 2017, Data-driven multi-task sparse dictionary learning for noise attenuation of 3D seismic data, Geophysics, 82, V385-V396.
% Siahsar, M. A. N., V. Abolghasemi, and Y. Chen, 2017, Simultaneous denoising and interpolation of 2D seismic data using data-driven non-negative dictionary learning, Signal Processing, 141, 309-321.
% Chen, Y., M. Zhang, M. Bai, and W. Chen, 2019, Improving the signal-to-noise ratio of seismological datasets by unsupervised machine learning, Seismological Research Letters, 90, 1552-1564.
% Chen, Y., S. Zu, W. Chen, M. Zhang, and Z. Guan, 2019, Learning the blending spikes using sparse dictionaries, Geophysical Journal International, 218, 1379?1397. 
% Wang, H., Q. Zhang, G. Zhang, J. Fang, and Y. Chen, 2020, Self-training and learning the waveform features of microseismic data using an adaptive dictionary, Geophysics, 85, KS51?KS61.
% Zu, S., H. Zhou, R. Wu, M. Jiang, and Y. Chen, 2019, Dictionary learning based on dip patch selection training for random noise attenuation, Geophysics, 84, V169?V183.
% Zu, S., H. Zhou, R. Wu, and Y. Chen, 2019, Hybrid-sparsity constrained dictionary learning for iterative deblending of extremely noisy simultaneous-source data, IEEE Transactions on Geoscience and Remote Sensing, 57, 2249-2262.
% etc. 

if nargin==2
   mode=1; 
end

g = double(g); % in case of data format is unit8,12,16
f = double(f);

if ndims(f)<ndims(g)
    error ('Dimesion of two images don''t match!');
end

if mode ==1
s = size(g,3);
if s==1 % single channel    
    psnr = 20.*log10(norm(g(:,:),'fro')/norm(g(:,:)-f(:,:),'fro'));   
else % multi-channel
    
    psnr = zeros(s,1);
    for i = 1:s
        psnr(i) = 20.*log10(norm(g(:,:,i),'fro')/norm(g(:,:,i)-f(:,:,i),'fro'));
    end
end
else
    [n1,n2,n3]=size(g);
     psnr = 20.*log10(norm(reshape(g,n1,n2*n3),'fro')/norm(reshape(g-f,n1,n2*n3),'fro'));   
end
end
