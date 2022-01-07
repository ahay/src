 function FXEMD(infile,outfile,n1,n2,dt,lf,hf,N,verb)
% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Oct, 2014

% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API
    
%  Copyright (C) 2014 The University of Texas at Austin
%  Copyright (C) 2014 Yangkang Chen
%  
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%  
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  
%  REFERENCES:
%  Chen, Y. and J. Ma, 2014, Random noise attenuation by f-x empirical mode decomposition predictive filtering, Geophysics, 79, V81-V91.
%  Chen, Y., C. Zhou, J. Yuan, and Z. Jin, 2014, Application of empirical mode decomposition to random noise attenuation of seismic data, Journal of Seismic Exploration, 23, 481-495.
%  Chen, Y., S. Gan, T. Liu, J. Yuan, Y Zhang, and Z. Jin, 2015, Random noise attenuation by a selective hybrid approach using f-x empirical mode decomposition, Journal of Geophysics and Engineering, 12, 12-25.
%  Chen, Y., G. Zhang, S. Gan, and C. Zhang, 2015, Enhancing seismic reflections using empirical mode decomposition in the flattened domain, Journal of Applied Geophysics, 119, 99-105.
%  Gan, S., S. Wang, Y. Chen, J. Chen, W. Zhong, and C. Zhang, 2016, Improved random noise attenuation using f − x empirical mode decomposition and local similarity, Applied Geophysics, 13, 2016, 127-134.
%  Chen, Y., 2016, Dip-separated structural filtering using seislet thresholding and adaptive empirical mode decomposition based dip filter, Geophysical Journal International, 206, 457-469.
%  Chen, W., J. Xie, S. Zu, S. Gan, and Y. Chen, 2017, Multiple reflections noise attenuation using adaptive randomized-order empirical mode decomposition, IEEE Geoscience and Remote Sensing Letters, 14, 18-22.
%  Chen, Y., Y. Zhou, W. Chen, S. Zu, W. Huang, and D. Zhang, 2017, Empirical low rank decomposition for seismic noise attenuation, IEEE Transactions on Geoscience and Remote Sensing, 2017, 55, 4696-4711.
%  Chen, Y. and S. Fomel, 2018, EMD-seislet transform, Geophysics, 83, A27-A32.


%%from Madagascar to Matlab
% create memory
din=zeros(n1,n2);
rsf_read(din,infile)

%% Main program goes here !
dout=fxemd(din,lf,hf,dt,N,verb);
%% from Matlab to Madagascar

rsf_create(outfile,size(dout)');
rsf_write(dout,outfile);



