t,a,p,s,g,h=var('t,a,p,s,g,h')
tzo(t,a,g)= sqrt(t^2 +a^2*(1-g^2))
xzo(a,g)=a*(1-g^2)
f(t,p,h)=(t + sqrt(t^2 + 4*h^2*p^2))/2
t2(t,p,g,h)=t*f(t,p,h) + (1-1/g^2)*h^2
tp(t,p,s,h)=t2(t,p,s,h)/sqrt(f(t2(t,p,s,h),p*f(t,p,h),h))
xp(t,p,g,h)=h^2*p/f(t,p,h) * (1 - f(t,p,h)^2/f(t2(t,p,g,h),p*f(t,p,h),h))
tco(t,a,g,h)=tp(tzo(t,a,g),a/tzo(t,a,g),g,h)
xco(t,a,g,h)=xzo(a,g) + xp(tzo(t,a,g),a/tzo(t,a,g),g,h)
def Push(t,g):
    return parametric_plot([xco(t,a,g,1),g*tco(t,a,g,1)], (a,-3.5*t,3.5*t))
p=Push(0.8,0.8)+Push(1.2,0.8)+Push(1.6,0.8)+Push(2.,0.8)+Push(2.4,0.8)+Push(2.8,0.8)+Push(3.2,0.8)+Push(3.6,0.8)
p.axes_range(xmin=-4.75,xmax=4.75,ymin=7,ymax=0)
p.save(frame=True,axes=False,axes_labels=["Midpoint (km)","Pseudo-depth (km)"],fontsize=12,aspect_ratio=9.5/(7*golden_ratio),filename='junk_sage.pdf')