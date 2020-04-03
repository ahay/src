function [D,G]=yc_sgk(X,param)
%yc_sgk: SGK algorithm
% BY Yangkang Chen
% Jan, 2020
%
% INPUT
% X:     input training samples
% param: parameter struct
%   param.mode=1;   %1: sparsity; 0: error
%   param.niter=10; %number of SGK iterations to perform; default: 10
%   param.D=DCT;    %initial D
%   param.T=3;      %sparsity level
%
% OUTPUT
% D:    learned dictionary
% G:    sparse coefficients
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



T=param.T;%T=1; %requred by SGK
niter=param.niter;
mode=param.mode;
if isfield(param,'K')
    K=param.K;
else
    K=size(D,2);%dictionary size: number of atoms
end

D=param.D(:,1:K);

for iter=1:niter
    
    if mode==1
        G=ompN(D,X,1);
    else
        %error defined sparse coding
    end
    
    for ik=1:K %KSVD iteration, K times SVD
        inds=find(G(ik,:)~=0);
        if ~isempty(inds)
            D(:,ik)=sum(X(:,inds),2);
            D(:,ik)=D(:,ik)/norm(D(:,ik));
        end
    end
end

%% extra step
G=ompN(D,X,T);

return

function [ G ] = ompN( D, X, K )
%multi-column sparse coding
% BY Yangkang Chen
% Jan, 2020
[n1,n2]=size(X);
[n1,n3]=size(D);
G=zeros(n3,n2);
if K==1
    for i2=1:n2
        G(:,i2)=omp_e(D,X(:,i2));
    end
else
    for i2=1:n2
        G(:,i2)=yc_omp0(D,X(:,i2),K);
    end
end
return

function [ g ] = omp_e( D, x )
% yc_omp0: Most basic orthogonal matching pursuit for sparse coding
% this simple tutorial code use the sparseity-constrained sparse coding
% model
%
% two version of sparse coding problems:
% 1) Sparseity-constrained sparse coding problem
%   gamma = min |x-Dg|_2^2  s.t. |g|_0 <= K
% 2) Error-constrained sparse coding problem
%   gamma = min |g|_0       s.t. |x-Dg|_2^2 <=eps
%
% Author: Yangkang Chen
% Oct 25, 2016
[n1,n2]=size(D);
g=zeros(n2,1);

max=0;                      %initialize max=0
for i2=1:n2                 %loop over all atoms (greedy algorithm)
    dtr=abs(sum(D(:,i2).*x));
    if max<dtr
        max=dtr;
        k=i2;
    end
end

g(k)=sum(D(:,k).*x)/sum(D(:,k).*D(:,k));

return


function [ G2 ] = yc_spray( G, n )
%multi-column sparse coding
% BY Yangkang Chen
% Jan, 2020
% G: row or column vector
% axis: 1 or 2
% n:  size

[n1,n2]=size(G);
if n1==1 %row vector, axis=1
    G2=zeros(n,n2);
    for i1=1:n
        G2(i1,:)=G;
    end
else %column vector, axis=2
    G2=zeros(n1,n);
    for i2=1:n
        G2(:,i2)=G;
    end
end
return

function [ g ] = yc_omp0( D, x, K )
% yc_omp0: Most basic orthogonal matching pursuit for sparse coding
% this simple tutorial code use the sparseity-constrained sparse coding
% model
% 
% two version of sparse coding problems:
% 1) Sparseity-constrained sparse coding problem
%   gamma = min |x-Dg|_2^2  s.t. |g|_0 <= K 
% 2) Error-constrained sparse coding problem
%   gamma = min |g|_0       s.t. |x-Dg|_2^2 <=eps
% 
% Author: Yangkang Chen
% Oct 25, 2016
[n1,n2]=size(D);
I=[];
r=x;
g=zeros(n2,1);
for ik=1:K
    k=[];
    max=0;                      %initialize max=0
    for i2=1:n2                 %loop over all atoms (greedy algorithm)       
        if sum(find(I==i2))==0  %search among the other atoms
            dtr=abs(sum(D(:,i2).*r));
            if max<dtr
                max=dtr;
                k=i2;
            end        
        end
    end
    I=[I,k];
    g(I)=inv(D(:,I)'*D(:,I))*D(:,I)'*x; %g_I = D_I^{+}x, D_I^TD_I is guaranteed to be non-singular 
    r=x-D(:,I)*g(I);
end

return
