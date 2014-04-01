t,p,s,g,h,z=var('t,p,s,g,h,z')
f(t,p,h)=(t + sqrt(t^2 + 4*h^2*p^2))/2
t2(t,p,g,h)=t*f(t,p,h) + (1-1/g^2)*h^2
t0(t,p,g,h)=f(t2(t,p,g,h),p*f(t,p,h),h)
x0(t,p,h)=h^2*p/f(t,p,h)
xr(z,t,p,g,h)=x0(t,p,h) - t0(t,p,g,h)/(p*f(t,p,h)) * (1 - z/t2(t,p,g,h))
tp(t,p,s,h)=t2(t,p,s,h)/sqrt(f(t2(t,p,s,h),p*f(t,p,h),h))
xp(t,p,g,h)=h^2*p/f(t,p,h) * (1 - f(t,p,h)^2/f(t2(t,p,g,h),p*f(t,p,h),h))
def smile(t,g,h):
    return parametric_plot([x0(t,p,h),sqrt(t2(t,p,g,h))],(p,-2.1,2.1))
sm1=smile(1,1.2,1)
sm2=smile(1,1,1)
sm3=smile(1,0.8,1)
def triangle(t,g,h):
    return parametric_plot([xp(t,p,g,h),tp(t,p,g,h)],(p,-10,10),thickness=3,color='green')
tr1=triangle(1,1.2,1)
tr2=triangle(1,0.8,1)
def ray(p,t,g,h):
    return parametric_plot([xr(z,t,p,g,h),sqrt(z)],(z,tp(t,p,g,h)^2,t2(t,p,g,h)),linestyle='--',color='green') 
ra1=ray(-2,1,1.2,1)+ray(-1.75,1,1.2,1)+ray(-1.5,1,1.2,1)+ray(-1.25,1,1.2,1)+ray(-1,1,1.2,1)+ray(-0.75,1,1.2,1)+ray(-0.5,1,1.2,1)+ray(-0.25,1,1.2,1)+ray(2,1,1.2,1)+ray(1.75,1,1.2,1)+ray(1.5,1,1.2,1)+ray(1.25,1,1.2,1)+ray(1,1,1.2,1)+ray(0.75,1,1.2,1)+ray(0.5,1,1.2,1)+ray(0.25,1,1.2,1)
ra2=ray(-2,1,1,1)+ray(-1.75,1,1,1)+ray(-1.5,1,1,1)+ray(-1.25,1,1,1)+ray(-1,1,1,1)+ray(-0.75,1,1,1)+ray(-0.5,1,1,1)+ray(-0.25,1,1,1)+ray(2,1,1,1)+ray(1.75,1,1,1)+ray(1.5,1,1,1)+ray(1.25,1,1,1)+ray(1,1,1,1)+ray(0.75,1,1,1)+ray(0.5,1,1,1)+ray(0.25,1,1,1)
ra3=ray(-2,1,0.8,1)+ray(-1.75,1,0.8,1)+ray(-1.5,1,0.8,1)+ray(-1.25,1,0.8,1)+ray(-1,1,0.8,1)+ray(-0.75,1,0.8,1)+ray(-0.5,1,0.8,1)+ray(-0.25,1,0.8,1)+ray(2,1,0.8,1)+ray(1.75,1,0.8,1)+ray(1.5,1,0.8,1)+ray(1.25,1,0.8,1)+ray(1,1,0.8,1)+ray(0.75,1,0.8,1)+ray(0.5,1,0.8,1)+ray(0.25,1,0.8,1)
p1=sm1+tr1+ra1
p2=sm2+ra2
p3=sm3+tr2+ra3
p1.axes_range(xmin=-0.85,xmax=0.85,ymin=1.9,ymax=0.95)
p2.axes_range(xmin=-0.85,xmax=0.85,ymin=1.7,ymax=0.75)
p3.axes_range(xmin=-0.85,xmax=0.85,ymin=1.5,ymax=0.55)
p=graphics_array([[p1],[p2],[p3]])
p.save(frame=True,axes=False,axes_labels=['midpoint (km)','time (s)'],figsize=[6,12],filename='junk_sage.pdf')