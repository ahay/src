s,x,y=var('s,x,y')
spline(s,x,y) = 1/10*(114 - 73*s) + 8/15*(-11 + 8*s)*(cos(x) + cos(y)) - 7/15*(-1 + s)*(cos(2*x) + cos(2*y)) + 28/45*(-3 + 2*s)*(cos(x - y) + cos(x + y)) + 1/180*(-6 + 5*s)*(cos(2*(x - y)) + cos(2*(x + y))) - 4/45*(-9 + 8*s)*(cos(x-2*y)+cos(2*x-y)+cos(2*x+y)+cos(x+2*y))
def cp(s):
    return contour_plot(spline(s,x,y),(x,-pi,pi),(y,-pi,pi),contours=range(0,17,2))+text('tension=%g' % s,(0,0),color='white')
p = graphics_array([[cp(s) for s in (0,0.3)],[cp(s) for s in (0.7,1)]])
p.save(filename='junk_sage.pdf',frame=True,axes=False)