function Real(data,fxemdpf)
% Requirements: RSF (http://rsf.sourceforge.net/) with Matlab API

% parameters definition
flow=5;
fhigh=245;
dt=0.001;
N=3;
lf=10;
mu=0.01;
% twin=500;
% xwin=500;
twin=700;
xwin=266;
mode1=1;
mode2=2;
mode3=3;

% allocate memory
d = zeros(700,266);

% from Madagascar to Matlab
rsf_read(d,data);

% f-x emdpf denoise
d3= denoise_win(d,flow,fhigh,mode3,dt,N,lf,mu,twin,xwin);

% output 
rsf_create(fxemdpf,size(d3)');
rsf_write(d3,fxemdpf);
