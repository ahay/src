function psnr = yc_snr(g,f,mode)
% Author: Yangkang Chen
% g: ground truth image
% f: noisy/restored image
% mode:1->2D SNR, 2->3D SNR

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
