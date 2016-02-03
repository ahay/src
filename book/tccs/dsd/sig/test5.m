%% test5.m
% test2
dslet=zeros(200,256); % noisy
dslet1=zeros(200,256);% thr
rsf_read(dslet,'datan-slet.rsf');
rsf_read(dslet1,'datan-sletthr.rsf');
un=dslet;

% figure;imagesc([u,un]);
    lambda=0.5;
    niter=30;
    lvl=2;
    htype='spline';
    thr=0.02;

u1=ddtf_denoise2d(un, lambda, niter, lvl, htype, thr);

figure;imagesc([un,u1,un-u1],[0,1]);

rsf_write(u1,'datan-sletddtfthr.rsf');
