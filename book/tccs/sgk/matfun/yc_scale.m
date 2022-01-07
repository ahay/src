function [ D1 ] = yc_scale(D,N,dscale)
%yc_scale: Scale the data up to the Nth dimension = sfscale axis=N
%  IN   D:   	intput data
%       N:      number of dimension for scaling
%               default: N=2
% 		dscale:  Scale by this factor
%       (does not include the rscale and pclip functions (not convenient actually))
%
%  OUT   D1:  	output data
%
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
%  Modified by Yangkang Chen on Jan, 2020
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
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


if nargin==0
    error('Input data must be provided!');
end

if nargin==1
    N=2;
    dscale=1.0;
end

if nargin==2
    dscale=1.0;
end


[n1,n2,n3]=size(D);
D1=D;
switch N
    case 1
        for i3=1:n3
            for i2=1:n2
                D1(:,i2,i3)=D1(:,i2,i3)/max(abs(D1(:,i2,i3)));
            end
        end
    case 2
        for i3=1:n3
            D1(:,:,i3)=D1(:,:,i3)/max(max(abs(D1(:,:,i3))));
        end
    case 3
        D1=D1/max(max(max(abs(D1))));
        
    case 0
        
        D1=D1*dscale;
        
    otherwise
        error('Invalid argument value N.')
end


return
