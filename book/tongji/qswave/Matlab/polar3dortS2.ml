% Note: We follow Yan & Sava (2009)'s work to calculate the polarization
%       vectors for qP, qSH and qSV modes

threshold=0.002; 
% "Standard" Orthorhombic Model (Schoenberg & Helbig, 1997)
%thoms[vp0,vs0,eps1,eps2,del1,del2,del3,gam1,gam2] %
thoms = [2437,1265,0.329,0.258,0.083,-0.078,-0.106,0.182,0.0455];
cc = thoms2stiff_ort(thoms);
bipolar3dtest(cc,threshold,'s2');
print -depsc junk_ml.eps;
