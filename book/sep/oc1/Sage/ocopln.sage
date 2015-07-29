x,y,h=var('x,y,h')
ymax(x,h)= 1/2*(x + sqrt(4*h^2 + x^2))		
def f2(h):
    return plot(sqrt(y^2-h^2), (y, h, ymax(5.1, h)),thickness=3)	
fronts2  = add(map(f2,range(6)))		
def r2(x):
    return plot(sqrt(y*x), (y, 0, ymax(x, 5)),linestyle='--',color='red') 
rays2 = add(map(r2,range(1,6)))
plot2  =  rays2+fronts2
plot2.axes_range(ymin=12,ymax=0)  
def f1(h): 
    return plot(sqrt(x^2 + 3*h^2), (x, ymax(0, h), ymax(5.1, h)),thickness=3) 
fronts1  = add(map(f1,range(6)))
def r1(y):
    return plot(sqrt(4*x^2 - 3*x*y), (x, 3/4*y, ymax(y, 5)),linestyle='--',color='red') 
rays1 = add(map(r1,range(6)))
plot1  =  rays1+fronts1
plot1.axes_range(ymin=12,ymax=0) 
plots=graphics_array([plot1,plot2])
plots.save(frame=True,axes_labels=['midpoint','time'],fontsize=14,figsize=[12,6],filename='junk_sage.pdf')