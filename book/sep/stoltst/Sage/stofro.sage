j,q,t=var('j,q,t')
x(j,q,t) = (j*cos(t))/sqrt(q)
z(j,q,t) = (j*sin(t))/q + j*(1/q - 1) + 1 
def g(q,a):
    return parametric_plot([x(a, q, t), z(a, q, t)], (t, 0, 2*pi))  
def gplot(q):
    p = g(q, 0.001) + g(q, 0.2) + g(q, 0.4) + g(q, 0.6) + g(q, 0.8)
    p.axes_range(xmin=-1.2,xmax=1.2,ymin=4,ymax=0)
    return p
p=graphics_array([gplot(1/2),gplot(3/2)])
p.save(frame=True,axes=False,filename='junk_sage.pdf')