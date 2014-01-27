a,y=var('a,y')
f(a,y)=(2*(-(1 + a^2*cos(y)^2 + tan(y)^2 + (a^4*cos(y)^4*sin(y)^2)/(1 + a^2*tan(y)^2))^(1/2) + ((1 + (-a + tan(y))^2)^(1/2) + (1 + (a + tan(y))^2)^(1/2))/2))/((1 + (-a + tan(y))^2)^(1/2) + (1 + (a + tan(y))^2)^(1/2))
def lo(a):
    return plot(f(a,y*pi/180),y,0,90,axes_labels=['angle','relative error'])
p=lo(1)
p.save(filename='junk_sage.pdf')