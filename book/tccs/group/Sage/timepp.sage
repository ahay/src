c11,c12,c13,c22,c23,c33,c44,c55,c66=var('c11,c12,c13,c22,c23,c33,c44,c55,c66')
M = matrix([
[c11, c12, c13, 0, 0, 0],
[c12, c22, c23, 0, 0, 0],
[c13, c23, c33, 0, 0, 0],
[0, 0, 0, c44, 0, 0],
[0, 0, 0, 0, c55, 0],
[0, 0, 0, 0, 0, c66]])

n1,n2,n3=var('n1,n2,n3') # components of the normal vector
nv=vector([n1,n2,n3]).column()

A = matrix([
[n1, 0, 0, 0, n3, n2],
[0, n2, 0, n3, 0, n1],
[0, 0, n3, n2, n1, 0]])
V=M.substitute(c12=c11-2*c66,c22=c11,c23=c13,c44=c55)
C=A*V*A.transpose()
e3=map(lambda x: x.full_simplify(), C.eigenvalues())

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

vgptrue(n1)=vgs(n1)
vgstrue(n1)=vgp(n1)

sn12(n1)=(vsgroup[0]^2/vgs).full_simplify()
pn12(n1)=(vpgroup[0]^2/vgp).full_simplify()

pn12true(n1) = sn12(n1)


# Muir-Dellinger
Q,N1=var('Q,N1')

ELp(n1)=(1/c11)*pn12(n1)+(1/c33)*(1-pn12(n1))
ELs(n1)=(1/c11)*sn12(n1)+(1/c33)*(1-sn12(n1))

qz=((2*c13 + c33)*c55 + c13^2)/(c11*c33 - c11*c55)
QZ = 1/qz

# Shifted hyperbola

S=var('S')
SHp(n1)=((1-S)*ELp(n1)+S*sqrt(ELp(n1)^2+2*(Q-1)*(1/c11)*(1/c33)*pn12*(1-pn12)/S))^-1
SHs(n1)=((1-S)*ELs(n1)+S*sqrt(ELs(n1)^2+2*(Q-1)*(1/c11)*(1/c33)*sn12*(1-sn12)/S))^-1

SHptrue(n1)=SHs(n1)

gsin(x)=(sqrt(pn12true(sin(x)))).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51) # sine of group angle
gcos(x)=sqrt(1-gsin(x)*gsin(x)) # cosine of group angle

TS(n1)= (1-sn12(n1))/c33 + (sn12(n1)*(-c55 + c33)*(sn12(n1)*c33*(-c55 + c33) + (1-sn12(n1))*(c13^2 + c55*(c33 + 2*c13))))/(c11*sn12(n1)*(c55 - c33)^2*c33 + (1-sn12(n1))*(c13^2 + c55*(c33 + 2*c13))^2)

SHplot = parametric_plot([gsin(x)/gcos(x),2000/gcos(x)*abs(1/sqrt(SHptrue(sin(x)))-1/sqrt(vgptrue(sin(x)))).subs(S=1/(2*(1+QZ)),Q=QZ).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)],(x,0,1.33))
ZUplot = parametric_plot([gsin(x)/gcos(x),2000/gcos(x)*abs(1/sqrt(SHptrue(sin(x)))-1/sqrt(vgptrue(sin(x)))).subs(S=1/2,Q=QZ).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)],(x,0,1.33),linestyle='--',color='green')
TSplot = parametric_plot([gsin(x)/gcos(x),2000/gcos(x)*abs(sqrt(TS(sin(x)))-1/sqrt(vgptrue(sin(x)))).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51)],(x,0,1.33),linestyle=':',color='red')

p=SHplot+ZUplot+TSplot
p.save(filename='junk_sage.pdf',frame=True,axes_labels=['offset (km)','absolute time error (ms)'],aspect_ratio=1/4,axes=False)
