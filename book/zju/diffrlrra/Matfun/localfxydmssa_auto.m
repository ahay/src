function [ DATA_out ] = localfxydmssa_auto( DATA_in, param )
% Local DRR with automatically selected ranks
%
%  Copyright (C) 2018 Yangkang Chen
% 
% DATA_in:      input data
% DATA_out:     output data
% param:        parameters set
% param.dt:     sampling interval in sec
% param.lf:   	min  freq. in the data in Hz
% param.hf:  	max  freq. in the data in Hz
% param.N  	
%
%  REFERENCES
%  Wang et al., 2020, Separation and imaging of seismic diffractions using a localized rank-reduction method with adaptively selected ranks, doi: 10.1190/geo2020-0215.1.
%  Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.
%  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
%  Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
%  Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
%  Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
% 
dt=param.dt;
N=param.N;
NN=param.NN;
flow=param.flow;
fhigh=param.fhigh;
verb=param.verb;
mode=param.mode;
f=0;

DATA_out=fxydmssa_auto(DATA_in,flow,fhigh,dt,N,NN,verb,mode,f,eps);

return


