o=var('o')
ph(o)=sqrt(2.0)*gamma(0.5-I*o*0.5)/gamma(-I*o*0.5)
ap(o)=(-I*o)^0.5
p1=plot(abs(ph(o)),(o,0,10))+plot(abs(ap(o)),(o,0,10),linestyle='--',color='green')+text('$|F|$',(5,3),fontsize=16,color='black')
myarg(x)=arctan(imag(x)/real(x))/pi
p2=plot(myarg(ph(o)),(o,0.001,10),ymin=-0.5,ymax=-0.2)+plot(myarg(ap(o)),(o,0.001,10),linestyle='--',color='green',ymin=-0.5,ymax=-0.2)+text(r'$\arg(F)/\pi$',(5,-0.22),fontsize=16,color='black')
p=graphics_array([p1,p2])
p.save(axes_labels=['$\omega$',''],fontsize=16,figsize=[12,6],filename='junk_sage.pdf')