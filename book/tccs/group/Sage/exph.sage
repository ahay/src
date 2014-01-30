c11,c12,c13,c22,c23,c33,c44,c55,c66=var('c11,c12,c13,c22,c23,c33,c44,c55,c66')
M = matrix([
[c11, c12, c13, 0, 0, 0],
[c12, c22, c23, 0, 0, 0],
[c13, c23, c33, 0, 0, 0],
[0, 0, 0, c44, 0, 0],
[0, 0, 0, 0, c55, 0],
[0, 0, 0, 0, 0, c66]])

n1,n2,n3=var('n1,n2,n3') # components of the normal vector
A = matrix([
[n1, 0, 0, 0, n3, n2],
[0, n2, 0, n3, 0, n1],
[0, 0, n3, n2, n1, 0]])

V=M.substitute(c12=c11-2*c66,c22=c11,c23=c13,c44=c55)
C=A*V*A.transpose()
e3=map(lambda x: x.full_simplify(), C.eigenvalues())

vp(n1)=e3[1].substitute(n2=0,n3=sqrt(1-n1^2))
vsv(n1)=e3[2].substitute(n2=0,n3=sqrt(1-n1^2))

# Sage bug : Always pick the negative root under the square root so we swictch the expression of qP and qSV
vptrue=vsv
vsvtrue=vp

ppx=(sqrt(vptrue(sin(x)))*sin(x)).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)
ppz=(sqrt(vptrue(sin(x)))*cos(x)).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)
pp=parametric_plot([ppx,ppz],(x,0,2*pi))

psvx=(sqrt(vsvtrue(sin(x)))*sin(x)).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)
psvz=(sqrt(vsvtrue(sin(x)))*cos(x)).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)
psv=parametric_plot([psvx,psvz],(x,0,2*pi),color='green')
p=pp+psv

p.save(filename='junk_sage.pdf',frame=True,axes_labels=['horizontal velocity (km/s)','vertical velocity (km/s)'])
