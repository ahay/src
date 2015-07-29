x1,x2,z1,z2=var('x1,x2,z1,z2')
R1,R2,R3=var('R1,R2,R3')
x,z,x3,z3=var('x,z,x3,z3')

eq1=(x-x2)^2+(z-z2)^2==R2^2
eq2=(x2-x)*(x3-x1)+(z2-z)*(z3-z1)==R2*(R3-R1)

def c(x,z,r,a):
    return circle((x,z),r)+circle((x,z),0.2,fill=True)+text(a,(x+0.5,z+0.5),fontsize=14)
def cr(x,z,r):
    return circle((x,z),r,thickness=2)

ln(x,x1,x2,z1,z2)=x*(z1-z2)/(x1-x2) + (z2*x1 - z1*x2)/(x1-x2)
cline = plot(ln(x,1, 6, 1, 0),(x,-1,10), linestyle='--') 

zs=solve(eq2,z)
sline=z.subs(zs[0]).subs(R1=1,R2=2,R3=3,x1=1,x2=3,x3=6,z1=1,z2=-1,z3=0)
eline = plot(sline,(x,1,3),thickness=2)
sline=z.subs(zs[0]).subs(R1=1,R2=0,R3=3,x1=1,x2=3,x3=6,z1=1,z2=-1,z3=0)
aline = plot(sline,(x,2,4))

xzs=solve([eq1,eq2],(x,z))
spoint = map(lambda r: r.subs(R1=1,R2=2,R3=3,x1=1,x2=3,x3=6,z1=1,z2=-1,z3=0),(x.subs(xzs[1][0]),z.subs(xzs[1][1])))
cpoint = circle(spoint,0.2,fill=True)+text('D',(spoint[0]-0.5,spoint[1]-0.5),fontsize=14)

p=c(1,1,1,'A')+c(3,-1,2,'B')+cr(3,-1,2)+c(6,0,3,'C')+eline+aline+cpoint+cline
p.save(axes=False,filename='junk_sage.pdf')