% Note: We follow Yan & Sava (2009)'s work to calculate the polarization
%       vectors for qP, qSH and qSV modes with some modification


threshold=0.002; 
% "Standard" Orthorhombic Model (Schoenberg & Helbig, 1997)
%thoms[vp0,vs0,eps1,eps2,del1,del2,del3,gam1,gam2] %
%thoms = [2437,1265,0.329,0.258,0.083,-0.078,-0.106,0.182,0.0455];
%cc = thoms2stiff_ort(thoms);
%bipolar3dtest(cc,threshold,'sh');
%print -depsc polar3dortSH.eps;

% VTI is a special case of ORT
% VTI model: Mesaverde shale (Thomsen, 1986)
%thoms = [3749,2621,0.225,0.078,0.100];
%cc = thoms2stiff_vti(thoms);
%bipolar3dtest(cc,threshold,'sh');
%print -depsc polar3dvtiSH.eps;

% HTI model: Isotropic background containing a vertical fratcure system
% with normal and tangential weaknesses 0.1 and 0.2 (Shang et al., 2015)
thoms = [3267,1697,-0.0447,-0.14,-0.100];
cc = thoms2stiff_hti(thoms);
bipolar3dtest_hti(cc,threshold,'sh');
print -depsc polar3dhtiSH.eps;

% ISO model: 
%thoms = [4721,2890];
%cc = thoms2stiff_iso(thoms);
%bipolar3dtest(cc,threshold,'sh');
%print -depsc polar3disoSH.eps;
