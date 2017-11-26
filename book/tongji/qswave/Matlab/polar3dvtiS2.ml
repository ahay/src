threshold=0.002; 

% VTI is a special case of ORT
% VTI model: Mesaverde shale (Thomsen, 1986)
thoms = [3749,2621,0.225,0.078,0.100];
cc = thoms2stiff_vti(thoms);
bipolar3dtest(cc,threshold,'s2');
print -depsc junk_ml.eps;
