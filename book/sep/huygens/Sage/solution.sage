x1,x2,z1,z2=var('x1,x2,z1,z2')
R1,R2,R3=var('R1,R2,R3')
x,z,x3,z3=var('x,z,x3,z3')

eq1=(x-x2)^2+(z-z2)^2==R2^2
eq2=(x2-x)*(x3-x1)+(z2-z)*(z3-z1)==R2*(R3-R1)

# definitions for various graphical objects
# circle with disk
def c(x,z,r,a):
    return circle((x,z),r)+circle((x,z),0.2,fill=True)+text(a,(x+0.5,z+0.5),fontsize=14)
# point
def p(x,z,a):
    return circle((x,z),0.2,fill=True)+text(a,(x+0.5,z+0.5),fontsize=14) 
def q(x,z,a):
     return circle((x,z),0.2,fill=True)+text(a,(x-0.5,z-0.5),fontsize=14)
# line
ln(x,x1,z1,x2,z2)=x*(z1-z2)/(x1-x2) + (z2*x1 - z1*x2)/(x1-x2)

# line connecting A with C
cline = plot(ln(x,1, 1, 6, 0),(x,1,6),linestyle='--') 

# solve for z=f(x) - draw eq2

zs=solve(eq2,z)
sline=z.subs(zs[0]).subs(R1=1,R2=2,R3=3,x1=1,x2=3,x3=6,z1=1,z2=-1,z3=0)
eline = plot(sline,(x,1,3.5)) 

#solve system

xzs=solve([eq1,eq2],(x,z))
spoint1 = map(lambda r: r.subs(R1=1,R2=2,R3=3,x1=1,x2=3,x3=6,z1=1,z2=-1,z3=0),(x.subs(xzs[0][0]),z.subs(xzs[0][1])))
spoint2 = map(lambda r: r.subs(R1=1,R2=2,R3=3,x1=1,x2=3,x3=6,z1=1,z2=-1,z3=0),(x.subs(xzs[1][0]),z.subs(xzs[1][1])))

# plot D and E - solutions of the system
dpoint=q(spoint1[0],spoint1[1],"D")
epoint=q(spoint2[0],spoint2[1],"E")

#plot OD and OE
odline = plot(ln(x1, 5, 4, spoint1[0], spoint1[1]),(x1,spoint1[0],5),thickness=2)
oeline = plot(ln(x1, 5, 4, spoint2[0], spoint2[1]),(x1,spoint2[0],5),thickness=2)

#tangents
#{z1 == a x1 + b, 
# z2 == a x2 + b, 
#(z1 - cz1)^2 == r1^2 - (x1-cx1)^2,
#(z2 - cz2)^2 == r2^2 - (x2-cx2)^2,
#(x1 - cx1) (x1-x2) + (z1 - cz1) (z1 - z2) == 0,
#(x2 - cx2) (x1-x2) + (z2 - cz2) (z1 - z2) == 0
#};
#%/.{cx1->1,cx2->6,cz1->1,cz2->0,r1->1,r2->3};
#stang = NSolve[%];
#ctang = Plot[(a x + b)/.stang[[2]], {x,-1,8}, \
#PlotStyle->Dashing[{0.05,0.05}]];
#a x+c == z/.stang[[2]];
#%/.spoint[[1]];
#Solve[%,c];
#tang = Plot[(a x + c)/.stang[[2]]/.%[[1]], {x,-1,8}];

opoint=p(5,4,"O")
apoint=c(1,1,1,"A")
bpoint=c(3,-1,2,"B")
cpoint=c(6,0,3,"C")

picture = opoint + apoint + bpoint + cpoint + dpoint + epoint + cline + eline + odline + oeline
#ctang, tang

picture.save(axes=False,filename='junk_sage.pdf')