% Author      : Yangkang Chen
%               Texas Consortium of Computational Seismology
%               Jackson School of Geosciences
%               The University of Texas at Austin
%         
% Date        : Apr, 2014 

clear;clc;close all;
%% generate hyper bolic events (no conflicts)
  dt = 2./1000;
  tmax = 1.0;
  h = [-500:20:1000];
  tau = [0.1,.65,0.85];
  v = [1500,2400,2300];
  amp = [1., -1.,1];
  f0 = 20;
  snr = 200;
  L = 9;
  seed=2013;
  
hyper=hyperbolic_events(dt,f0,tmax,h,tau,v,amp,snr,L,seed);
hyper=hyper/max(max(hyper));

fid=fopen('../hyper/hyper.bin','w');
fwrite(fid,hyper,'float');
%figure;imagesc(d);colormap(seismic);

%% generate complex linear events (with conflicts)
linear = linear_events(0.004,40,2,[0:10:10*79],[0.25,0.5,1.0,1.35],[0,0.0006,-0.0006,0.0006],[1,1,1,1],200,2,201314);
linear=linear/max(max(linear));
fid=fopen('../complex/complex.bin','w');
fwrite(fid,linear,'float');
%figure;imagesc(linear);colormap(seismic);


















