t,p,s,h,v=var('t,p,s,h,v')
f(t,p,h)=(t + sqrt(t^2 + 4*h^2*p^2))/2
t2(t,p,s,h)=t*f(t,p,h) + s*h^2
tp(x,p,s,h)=t2(p*x,p,s,h)/sqrt(f(t2(p*x,p,s,h),p*f(p*x,p,h),h))
xp(x,p,s,h)=x + h^2*p/f(p*x,p,h) * (1 - f(p*x,p,h)^2/f(t2(p*x,p,s,h),p*f(p*x,p,h),h))
def fp(v):
    return parametric_plot([xp(x,0.5,1-1/v^2,2),v*tp(x,0.5,1-1/v^2,2)],(x,0,2.1),thickness=3)
fps=add(map(fp,[0.8,0.9,1.0,1.1,1.2]))
def rp(x):
    return parametric_plot([xp(x,0.5,1-1/v^2,2),v*tp(x,0.5,1-1/v^2,2)],(v,0.8,1.2),linestyle='--',color='red')
rps=add(map(rp,[0.0,0.5,1.0,1.5,2.0]))
p1=fps+rps
p1.axes_range(xmin=0,xmax=2.5,ymin=2,ymax=0)

d(x)=sqrt(1 + x^2)
td(x,s,h)=t2(d(x),x/d(x),s,h)/sqrt(f(t2(d(x),x/d(x),s,h),x/d(x)*f(d(x),x/d(x),h),h))
xd(x,s,h)=x + h^2*x/(d(x)*f(d(x),x/d(x),h)) * (1 - f(d(x),x/d(x),h)^2/f(t2(d(x),x/d(x),s,h),x/d(x)*f(d(x),x/d(x),h),h))
def fd(v):
    return parametric_plot([xd(x,1-1/v^2,2),v*td(x,1-1/v^2,2)],(x,-2.1,2.1),thickness=3)
fds=add(map(fd,[0.90,0.95,1.00,1.05,1.10]))
def rd(x):
    return parametric_plot([xd(x,1-1/v^2,2),v*td(x,1-1/v^2,2)],(v,0.9,1.1),linestyle='--',color='red')
rds=add(map(rd,[-1.75,-1.25,-0.75,-0.25,0.25,0.75,1.25,1.75]))
p2=fds+rds
p2.axes_range(xmin=-2.5,xmax=2.5,ymin=3,ymax=0)

p=graphics_array([p2,p1])
p.save(fontsize=14,figsize=[12,6],filename='junk_sage.pdf')