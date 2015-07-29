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

nv=vector([n1,n2,n3]).column()
vsphase=sqrt(e3[2])
vpphase=sqrt(e3[1])

vsgroup=vsphase*nv+(identity_matrix(3)-nv*nv.transpose())*vector([diff(vsphase,n1),diff(vsphase,n2),diff(vsphase,n3)]).column()
vsgroup=[vsgroup[0][0].full_simplify(),vsgroup[1][0].full_simplify(),vsgroup[2][0].full_simplify()]
vsgroup=map(lambda x: x.substitute(n2=0,n3=sqrt(1-n1^2)).full_simplify(), vsgroup)

vpgroup=vpphase*nv+(identity_matrix(3)-nv*nv.transpose())*vector([diff(vpphase,n1),diff(vpphase,n2),diff(vpphase,n3)]).column()
vpgroup=[vpgroup[0][0].full_simplify(),vpgroup[1][0].full_simplify(),vpgroup[2][0].full_simplify()]
vpgroup=map(lambda x: x.substitute(n2=0,n3=sqrt(1-n1^2)).full_simplify(), vpgroup)

vgs(n1)=(vsgroup[0]^2+vsgroup[1]^2+vsgroup[2]^2).full_simplify()
vgp(n1)=(vpgroup[0]^2+vpgroup[1]^2+vpgroup[2]^2).full_simplify()

# Sage bug : Always pick the negative root under the square root so we swictch the expression of qP and qSV
vgptrue=vgs
vgstrue=vgp

sn12(n1)=(vsgroup[0]^2/vgs).full_simplify()
pn12(n1)=(vpgroup[0]^2/vgp).full_simplify()

pn12true=sn12
sn12true=pn12

vsx=(sqrt(vgstrue(sin(x)))*sign(sin(x))*sqrt(sn12true(sin(x)))).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)
vsz=(sqrt(vgstrue(sin(x)))*sign(cos(x))*sqrt(1-sn12true(sin(x)))).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)
pgs=parametric_plot([vsx,vsz],(x,0,2*pi),color='green')

vpx=(sqrt(vgptrue(sin(x)))*sign(sin(x))*sqrt(pn12true(sin(x)))).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)
vpz=(sqrt(vgptrue(sin(x)))*sign(cos(x))*sqrt(1-pn12true(sin(x)))).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)
pgp=parametric_plot([vpx,vpz],(x,0,2*pi))

p=pgp+pgs

p.save(filename='junk_sage.pdf',frame=True,axes_labels=['horizontal velocity (km/s)','vertical velocity (km/s)'])
