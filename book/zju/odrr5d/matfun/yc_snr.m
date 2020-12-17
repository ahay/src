function psnr = yc_snr(g,f,mode)
% Author: Yangkang Chen
% g: ground truth image
% f: noisy/restored image
% mode:1->2D SNR, 2->3D SNR
%
%  References:   
%
%  [1] Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  [2] Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  [3] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  [4] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  [5] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  [6] Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
%  [7] Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.

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
