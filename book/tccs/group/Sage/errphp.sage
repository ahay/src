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

# Muir-Dellinger approximation
q=var('q')
el(n1)=c11*n1^2+c33*(1-n1^2)
md(n1)=el(n1)+(q-1)*c11*c33*n1^2*(1-n1^2)/el(n1)
qz=((2*c13 + c33)*c55 + c13^2)/(c11*c33 - c11*c55)
mdplot=plot(100*abs(sqrt(md(sin(x*pi/180))/vptrue(sin(x*pi/180)))-1).substitute(q=qz).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51),x,0,90,linestyle='--')

# Shifted hyperbola approximation
s=var('s')
sh(n1)=(1-s)*el(n1)+s*sqrt(el(n1)^2+2*(q-1)*c11*c33*n1^2*(1-n1^2)/s)
shplot=plot(100*abs(sqrt(sh(sin(x*pi/180))/vptrue(sin(x*pi/180)))-1).substitute(q=qz).substitute(s=0.5).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51),x,0,90)

#weak anisotropy approximation
epsilon,delta=var('epsilon,delta')
epsilon=(c11-c33)/(2*c33)
delta=((c55+c13)^2-(c33-c55)^2)/(2*c33* (c33-c55))
th(n1)=c33*(1+2*epsilon*n1^4+2*delta*n1^2*(1-n1^2))
weakplot=plot(100*abs(sqrt(th(sin(x*pi/180))/vptrue(sin(x*pi/180)))-1).substitute(c11=14.47,c33=9.57,c55=2.28,c13=4.51),x,0,90,linestyle=':')

p=shplot+mdplot+weakplot
p.save(filename='junk_sage.pdf',axes_labels=['phase angle (degrees)','relative error (%)'],aspect_ratio=50,frame=True,axes=False)