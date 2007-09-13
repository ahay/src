function eigen(in,out)
% EIGEN Eigenvalues of a matrix

mat = zeros(rsf_dim(in));
rsf_read(mat,in);
rsf_write(eig(mat),out);
