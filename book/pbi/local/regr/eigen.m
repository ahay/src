function eigen(in,out)
% EIGEN Eigenvalues of a matrix

dims = rsf_dim(in); 
mat = zeros(dims');
rsf_read(mat,in);
ev = abs(eig(mat));
rsf_write(ev,out);
