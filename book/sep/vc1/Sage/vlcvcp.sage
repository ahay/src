t,p,h,s,g=var('t,p,h,s,g')
f(t,p,h)=(t + sqrt(t^2 + 4*h^2*p^2))/2
t2(t,p,g,h)=t*f(t,p,h) + (1-1/g^2)*h^2
tp(t,p,s,h)=t2(t,p,s,h)/sqrt(f(t2(t,p,s,h),p*f(t,p,h),h))
xp(t,p,g,h)=h^2*p/f(t,p,h) * (1 - f(t,p,h)^2/f(t2(t,p,g,h),p*f(t,p,h),h))
def triangle(t,g,h):
    return parametric_plot([xp(t,p,g,h),tp(t,p,g,h)],(p,-10,10),thickness=3,color='green')
plot1=triangle(0.5,1.333,1)+triangle(1,1.333,1)+triangle(1.5,1.333,1)
plot1.axes_range(xmin=-0.5,xmax=0.5,ymin=2,ymax=0)
plot2=triangle(0.75,0.833,1)+triangle(1.25,0.833,1)+triangle(1.75,0.833,1)
plot2.axes_range(xmin=-0.5,xmax=0.5,ymin=2,ymax=0)
plot=graphics_array([plot1,plot2])
plot.save(frame=True,axes=False,gridlines=[False,True],axes_labels=['midpoint (km)','time (s)'],filename='junk_sage.pdf')