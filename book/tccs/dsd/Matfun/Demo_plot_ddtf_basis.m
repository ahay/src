clear all
close all
clc
%% read image
% name = 'lena.tif';
name = 'cameraman.png';
% name = 'barbara.png';
% name = 'data1.png';
% name = 'data2.png';
f = double(imread(name));

if size(f,3)>1
    f = double(rgb2gray(uint8(f)));
end



%% initial framelet
% level = 3; % filter size
% type = 'haar';
% type = 'dct';
type = 'spline'; level = 2;
[initialH, hsize] = initialFilter(type, level); % generate 2D filters

displayH(initialH); % display the filters
title('Initial Basis','fontsize',16,'fontweight','bold');
print(gcf, '-depsc','-r300','basis1.eps');
%% data driven framelet
lambda = 10;
maxits = 50;
% maxits = 100;
learnedH = ddTF(f, initialH,hsize, lambda,maxits);

displayH(learnedH);
title('Final Basis','fontsize',16,'fontweight','bold');
print(gcf, '-depsc','-r300','basis2.eps');

tlt = sprintf('\\lambda = %d, maxits = %d', lambda, maxits);
%title(tlt);

%% zero boundary condition
% d1 = ddtf_dec_z(f, initialH);
% d = ddtf_dec_z(f, learnedH);
% 
% tau = 0;
% 
% if strcmp(type,'spline')
%     d1 = swt2(f,1,2);
%     d1_ = wavethresh(d1,tau,'soft');
%     r1 = iswt2(d1_,1,2);
% else
%     d1_ = ani_thresh(d1, 's', tau);
%     r1 = ddtf_rec_z(d1_, initialH);
% end
% 
% p1 = psnr(f, r1)
% 
% d_ = ani_thresh(d, 's', tau);
% r = ddtf_rec_z(d_, learnedH);
% 
% p2 = psnr(f, r)
%% circular boundary condition
% d1 = ddtf_dec_p(f, initialH);
% d = ddtf_dec_p(f, learnedH);
% 
% tau = 0;
% 
% if strcmp(type,'spline')
%     d1 = swt2(f,1,2);
%     d1_ = wavethresh(d1,tau,'soft');
%     r1 = iswt2(d1_,1,2);
% else
%     d1_ = ani_thresh(d1, 's', tau);
%     r1 = ddtf_rec_p(d1_, initialH);
% end
% 
% p1 = psnr(f, r1)
% 
% d_ = ani_thresh(d, 's', tau);
% r = ddtf_rec_p(d_, learnedH);
% 
% p2 = psnr(f, r)
