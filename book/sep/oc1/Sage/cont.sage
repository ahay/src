x,h=var('x,h')
t1(x,h)=sqrt(2.*(0.25 + h^2 - x^2 + sqrt(0.0625 - 0.5*h^2 + h^4 - 0.5*x^2 - 2.*h^2*x^2 + x^4)))
t2(x,h)=sqrt(2.*(0.25 + h^2 - x^2 - sqrt(0.0625 - 0.5*h^2 + h^4 - 0.5*x^2 - 2.*h^2*x^2 + x^4)))
def p1(h):
    return plot(t1(x,h),(x,-0.999999*abs(h-0.5),0.999999*abs(h-0.5)))
def p2(h):
    return plot(t2(x,h),(x,-0.999999*abs(h-0.5),0.999999*abs(h-0.5)))
pl1=p1(0.0)+p1(0.1)+p1(0.2)+p1(0.3)+p1(0.4)
pl2=p2(0.6)+p2(0.7)+p2(0.8)+p2(0.9)+p2(1.0)
pl=pl1+pl2
def scale_yaxis(numbers):
    return [-num for num in numbers]
pr1=p2(0.1)+p2(0.2)+p2(0.3)+p2(0.4)
pr2=p1(0.6)+p1(0.7)+p1(0.8)+p1(0.9)+p1(1.0)
pr=pl+pr1+pr2
p=graphics_array([pl,pr])
p.save(ymin=2.25,ymax=0,frame=True,aspect_ratio=1/2.25,axes_labels=['midpoint','time'],filename='junk_sage.pdf')