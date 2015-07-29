w,k,e=var('w,k,e')
fd(w,k) = (-27*I + w*(-12 + w*(-18*I + 5*w)) + (27*I + 12*w + w^3)*cos(k))/(-81*I + w*(24 + w*(-12*I + 5*w)) + (81*I + w*(-24 + w*(-6*I + w)))*cos(k))
F(e) = sqrt((1 + sqrt(1 + e^2))/(2*sqrt(1 + e^2)))*exp((1 - sqrt(1 + e^2))/2)
f(e) = 1/2*(1 - sqrt(1 + e^2) + log((1 + sqrt(1 + e^2))/2))
def asy(w,k):
    return -arg(F(2*k/w).n()/F(0).n()*exp(I*w*(f(2*k/w) - f(0)).n()))
def exe(w,k):
    return -arg(gamma(1 - (1 +  I*w)/2).n()*((k/2)^((1 +  I*w)/2)).n()*bessel_J(-(1 +  I*w)/2,k).n())
def exs(w):
    return [(x*pi/100,exe(w,x*pi/100)) for x in range(1,101)]
def ass(w):
    return [(x*pi/100,asy(w,x*pi/100)) for x in range(1,101)]
plot1 = line(exs(1),linestyle=':',color='black',legend_label='exact')+line(ass(1),linestyle='--',color='green',legend_label='asymptotic')+plot(-arg(fd(1,k)),(k,0,pi),legend_label='finite-difference')+text(r'$\Omega=1$',(0.5,0.8),axis_coords=True,fontsize=18) 
plot2 = line(exs(10),linestyle=':',color='black')+line(ass(10),linestyle='--',color='green')+plot(-arg(fd(10,k)),(k,0,pi))+text(r'$\Omega=10$',(0.5,0.8),axis_coords=True,fontsize=18) 
p = graphics_array([plot1,plot2])
p.save(frame=True,figsize=[12,6],axes_labels=['wavenumber (radians)',''],filename='junk_sage.pdf')