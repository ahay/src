t1=1
v=1
R = v*t1/2
h1=1
b = t1^2/(t1^2+4*h1^2/v^2)
xmax = R/sqrt(b)
f = pi/60
x,h=var('x,h')
x0(x)= x*(1-b)
dp(x,h)= -(-25 + 16*x^2 + sqrt((25 - 16*x^2)^2 + 1600*h^2*(x*cos(f) + 5*sqrt(1/4 - x^2/5)*sin(f))^2))/(40*(x*cos(f) + 5*sqrt(1/4 - x^2/5)*sin(f)))
dm(x,h)= -(-25 + 16*x^2 + sqrt((25 - 16*x^2)^2 + 1600*h^2*(x*cos(f) - 5*sqrt(1/4 - x^2/5)*sin(f))^2))/(40*(x*cos(f) - 5*sqrt(1/4 - x^2/5)*sin(f)))
tm(x,h)=t1*h/h1*sqrt((h1^2 - x0(x)^2)/(h^2 - dm(x,h)^2))
tp(x,h)=t1*h/h1*sqrt((h1^2 - x0(x)^2)/(h^2 - dp(x,h)^2))
ym(x,h)=x0(x)+dm(x,h)*cos(f)
yp(x,h)=x0(x)+dp(x,h)*cos(f)
def pp(h):
    return parametric_plot((yp(x,h),tp(x,h)),(x,-xmax,xmax),linestyle='--',color='black')+parametric_plot((ym(x,h),tm(x,h)),(x,-xmax,xmax),linestyle='--',color='black')
th(x,h)=t1*h*sqrt(2/(h^2 + h1^2 - x^2 + sqrt((h^2+h1^2-x^2)^2-4*h1^2*h^2)))
tl(x,h)=t1/h1*sqrt((h^2 + h1^2 - x^2 + sqrt((h^2+h1^2-x^2)^2-4*h1^2*h^2))/2)
def hp(h):
    return pp(h)+plot(th(x,h),(x,-(h-h1),h-h1))
def lp(h):
    return pp(h)+plot(tl(x,h),(x,-(h1-h),h1-h))
pl=hp(1.2)+hp(1.4)+hp(1.6)+lp(0.8)+lp(0.6)+lp(0.4)
pl.save(frame=True,axes=False,axes_labels=['midpoint','time'],filename='junk_sage.pdf')