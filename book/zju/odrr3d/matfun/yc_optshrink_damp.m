function [dout] = yc_optshrink_damp(din,N,K)
%yc_optshrink: Optimally shrink the singular values with damping
% BY Yangkang Chen, Sep, 2017
% Modified by Min Bai, May, 2018
% 
% INPUT
% din: input
% N:   rank
% 
% OUTPUT
% dout: output
%
% EXAMPLE
% test_odrr_denoise.m
%
% REFERENCE
%  Bai et al., 2020, Seismic signal enhancement based on the lowrank methods, Geophysical Prospecting, 68, 2783-2807.
% Chen et al., 2020, Five-dimensional seismic data reconstruction using the optimally damped rank-reduction method, Geophysical Journal International, 222, 1824-1845.
% Nadakuditi, R. R., 2013, Optshrink: An algorithm for improved low-rank signal matrix denoising by optimal, data-driven singular value shrinkage: IEEE Transactions on Information Theory, 60, 3002?3018.
% Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.
% Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
% Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
% Chen et al., 2017, Preserving the discontinuities in least-squares reverse time migration of simultaneous-source data, Geophysics, 82, S185-S196.
% Chen et al., 2019, Obtaining free USArray data by multi-dimensional seismic reconstruction, Nature Communications, 10:4434.

[n1,n2]=size(din);
[U,S,V]=svd(din);
s=diag(S);
if n1>=n2
   St=S(N+1:n2,N+1:n2); %make sure it is a square matrix
else
   St=S(N+1:n1,N+1:n1); %make sure it is a square matrix
end
weight=zeros(N,1);
for ir = 1 : N
    num=Ds(s(ir),St);
    den=DDs(s(ir),St);
    weight(ir) = -2*num/den*(1-(s(N+1)/s(ir))^K);
end

dout = U(:,1:N)*diag(weight,0)*V(:,1:N)';

end

function v = Ds(s,X)
% D transform of X
% INPUT
% s: singular value
% X: square singular value matrix
% OUTPUT
% v: output value
% 
% Reference
% eq 13 in Bai et al., 2020, GP. 

[n,m] = size(X);
if n~=m
    error('X must be a square matrix: n must be equal to m');
end
I=eye(m);
v=1./(m*m).*trace(s*inv(s^2*I-X.*X)).^2;

end

function v = DDs(s,X)
% derivative of D transform of X
% INPUT
% s: singular value
% X: square singular value matrix
% OUTPUT
% v: output value
% 
% Reference
% eq 14 in Bai et al., 2020, GP. 
[n,m] = size(X);
if n~=m
    error('X must be a square matrix: n must be equal to m');
end
I=eye(m);
v=2./(m*m).*trace(s*inv(s^2*I-X.*X)).*trace(inv(s^2*I-X.*X)-2*s*s*(inv(s*s*I-X.*X).^2));

end



