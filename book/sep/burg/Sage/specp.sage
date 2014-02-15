s,x,y=var('s,x,y')
spline(s,x,y) = 1/10*(114 - 73*s) + 8/15*(-11 + 8*s)*(cos(x) + cos(y)) - 7/15*(-1 + s)*(cos(2*x) + cos(2*y)) + 28/45*(-3 + 2*s)*(cos(x - y) + cos(x + y)) + 1/180*(-6 + 5*s)*(cos(2*(x - y)) + cos(2*(x + y))) - 4/45*(-9 + 8*s)*(cos(x-2*y)+cos(2*x-y)+cos(2*x+y)+cos(x+2*y))
exact(s,x,y) = (1-s)*(x^2+y^2)^2 + s*(x^2+y^2)
def pl(s):
    return plot(spline(s,x,0),(x,-pi,pi))+plot(exact(s,x,0),(x,-pi,pi),linestyle='--',color='black')+text('tension=%g' % s,(0,15),fontsize=12)
p = graphics_array([[pl(s) for s in (0,0.3)],[pl(s) for s in (0.7,1)]])
p.save(filename='junk_sage.pdf',ymax=16,frame=True,axes=False)