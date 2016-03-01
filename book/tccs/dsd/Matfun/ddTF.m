function [A] = ddTF(f, H,hsize, lambda,maxits)
% data driven framelet initialed by H
% avoid the transpose of large size matrix

[m,n] = size(f);
h_ = hsize - 1;

C = (m-h_)*(n-h_);
R = hsize^2;
mtx = zeros(R,C);
t = 1;
for i=1:m-h_
    for j=1:n-h_
        p = f(i:i+h_, j:j+h_);
        mtx(:,t) = p(:);
        t = t + 1;
    end
end

A = H;
G = mtx;

% maxits = 50;
for i = 1:maxits
    V = A'*G;
    V = wthresh(V, 'h', lambda);
    
    GV = G*V';
    [U, ~, X] = svd(GV);
    
    A = U*X';
    
end
