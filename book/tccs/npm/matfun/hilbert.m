function x = hilbert(xr,n)
%HILBERT  Discrete-time analytic signal via Hilbert transform.
%   X = HILBERT(Xr) computes the so-called discrete-time analytic signal
%   X = Xr + i*Xi such that Xi is the Hilbert transform of real vector Xr.
%   If the input Xr is complex, then only the real part is used: Xr=real(Xr).
%   If Xr is a matrix, then HILBERT operates along the columns of Xr.
%
%   HILBERT(Xr,N) computes the N-point Hilbert transform.  Xr is padded with 
%   zeros if it has less than N points, and truncated if it has more.  
%
%   For a discrete-time analytic signal X, the last half of fft(X) is zero, 
%   and the first (DC) and center (Nyquist) elements of fft(X) are purely real.
%
%   Example:
%     Xr = [1 2 3 4];
%     X = hilbert(Xr)
%   produces X=[1+1i 2-1i 3-1i 4+1i] such that Xi=imag(X)=[1 -1 -1 1] is the
%   Hilbert transform of Xr, and Xr=real(X)=[1 2 3 4].  Note that the last half
%   of fft(X)=[10 -4+4i -2 0] is zero (in this example, the last half is just
%   the last element).  Also note that the DC and Nyquist elements of fft(X)
%   (10 and -2) are purely real.
%
%   See also FFT, IFFT.

%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.10.4.3 $  $Date: 2007/12/14 15:05:05 $

%   References:
%     [1] Alan V. Oppenheim and Ronald W. Schafer, Discrete-Time
%     Signal Processing, 2nd ed., Prentice-Hall, Upper Saddle River, 
%     New Jersey, 1998.
%
%     [2] S. Lawrence Marple, Jr., Computing the discrete-time analytic 
%     signal via FFT, IEEE Transactions on Signal Processing, Vol. 47, 
%     No. 9, September 1999, pp.2600--2603.

if nargin<2, n=[]; end
if ~isreal(xr)
  warning(generatemsgid('Ignore'),'HILBERT ignores imaginary part of input.')
  xr = real(xr);
end
% Work along the first nonsingleton dimension
[xr,nshifts] = shiftdim(xr);
if isempty(n)
  n = size(xr,1);
end
x = fft(xr,n,1); % n-point FFT over columns.
h  = zeros(n,~isempty(x)); % nx1 for nonempty. 0x0 for empty.
if n > 0 && 2*fix(n/2) == n
  % even and nonempty
  h([1 n/2+1]) = 1;
  h(2:n/2) = 2;
elseif n>0
  % odd and nonempty
  h(1) = 1;
  h(2:(n+1)/2) = 2;
end
x = ifft(x.*h(:,ones(1,size(x,2))));

% Convert back to the original shape.
x = shiftdim(x,-nshifts);
