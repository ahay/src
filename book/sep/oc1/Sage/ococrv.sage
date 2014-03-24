y,x,h=var('y,x,h')
ymax(x,h) = (-1 + x^2 + sqrt(1 + (2 + 4*h^2)*x^2 + x^4))/(2*x)
ray(y,x) = sqrt(1 + x*y)
fro(y,h) = sqrt(1/2*(y^2  +  1  - h^2  + sqrt((y^2 + 1 - h^2)^2 + 4*h^2)))
def plotfro(h):
    return plot(fro(y, h), (y, -ymax(5.1, h), ymax(5.1, h)),thickness=3) 
fronts = add(map(plotfro,range(6)))
def plotrayp(x):
    return plot(ray(y, x), (y, -1/x, ymax(x, 5)),linestyle='--',color='red') 
def plotraym(x):
    return plot(ray(y, x), (y, ymax(x, 5), -1/x),linestyle='--',color='red')
rays = add(map(plotrayp,range(1,6))+map(plotraym,range(-5,0)))
plot1 = rays+fronts
plot1.axes_range(ymin=7,ymax=0) 
ray(y,x) = sqrt(1 - x*y)
def plotray1(x):
    return plot(ray(y, x), (y, -1, 1/x),linestyle='--',color='red')
def plotray2(x):
    return plot(ray(y, x), (y, 1/x, 1), linestyle='--',color='red')
rays2 = add(map(plotray1,[0.4,0.6,0.8,1.0])+map(plotray2,[-0.4,-0.6,-0.8,-1.0])) 
fro1(y,h) = sqrt(1/2*(1 + h^2 - y^2 + sqrt(-4*h^2 + (-1 - h^2 + y^2)^2)))
fro2(y,h) = sqrt(1/2*(1 + h^2 - y^2 - sqrt(-4*h^2 + (-1 - h^2 + y^2)^2)))
def plotfro1(h):
    return plot(fro1(y, h), (y, - (1 - h), 1 - h),thickness=3) 
def plotfro2(h):
    return plot(fro2(y, h), (y, - (h - 1), h - 1),thickness=3) 
fronts2=add(map(plotfro1,[0.8,0.6,0.4,0.2])+map(plotfro2,[1.2,1.4,1.6,1.8])) 
plot2=rays2+fronts2
plot2.axes_range(xmin=-1,xmax=1,ymin=1.4,ymax=0)
plots=graphics_array([plot1,plot2])
plots.save(frame=True,axes_labels=['midpoint','time'],fontsize=14,figsize=[12,6],filename='junk_sage.pdf')