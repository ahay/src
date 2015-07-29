def cir(x,z,r,a):
    return circle((x,z),r)+circle((x,z),0.1,fill=True)+ text(a,(x+0.5,z),fontsize=14)    
def cir2(x,z,a):
    return circle((x,z),0.1,fill=True)+text(a,(x+0.5,z),fontsize=14)  
x1,x2,z1,z2=var('x1,x2,z1,z2')
ln(x,x1,x2,z1,z2)=x*(z1-z2)/(x1-x2) + (z2*x1 - z1*x2)/(x1-x2)
abline = plot(ln(x,1, 4  , 1, -2),(x,1,4))
bcline = plot(ln(x,4, 4.5,-2,  2),(x,4,4.5))
obarrow = arrow((8,4),(4,-2),width=2)
p=cir2(8,4,'O')+cir(1,1,1,'A')+cir(4,-2,2,'B')+cir(4.5,2,3,'C')+abline+bcline+obarrow
p.save(axes=False,filename='junk_sage.pdf')