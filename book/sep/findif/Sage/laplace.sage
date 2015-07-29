x,y,g=var('x,y,g')
r(x,y,g)=2*(cos(sqrt(x^2+y^2)) - 1) - 2*(g*(cos(x) + cos(y) - 2) + (1-g)*(cos(x)*cos(y) - 1))
def p(g,a):
    return contour_plot(abs(r(x,y,g)),(x,0,pi),(y,0,pi),contours=srange(0.2,3.8,0.2),cmap='gray_r')+text(a,(0.1,0.1),axis_coords=True)
ps=graphics_array([p(1,'a'),p(0,'b'),p(0.5,'c'),p(0.65,'d')]) 
ps.save(filename='junk_sage.pdf',frame=True)