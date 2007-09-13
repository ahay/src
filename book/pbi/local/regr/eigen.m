function eigen(in,out)
% EIGEN Eigenvalues of a matrix

dims = rsf_dim(in); 
mat = zeros(dims');
rsf_read(mat,in);
rsf_write(eig(mat),out);
