def ccin(x):
    if abs(x) <= 1:
        return 3/2*abs(x)^3 - 5/2*abs(x)^2 + 1
    elif abs(x) > 1 and abs(x) <= 2: 
        return -1/2*abs(x)^3 + 5/2*abs(x)^2 - 4*abs(x) + 2
    else:
        return 0
w=var('w')
ccinf(w)= (2*integral(cos(w*x)*(3/2*x^3 - 5/2*x^2 + 1),(x,0,1)) + 2*integral(cos(w*x)*(-1/2*x^3 + 5/2*x^2 - 4*x + 2),(x,1,2))).full_simplify()
pccin=plot(ccin,(x,-2,2),thickness=3,aspect_ratio=4/golden_ratio)
wccin=plot(ccinf(w),(w,-2*pi,2*pi),ymin=-0.1,thickness=3,color='red',aspect_ratio=4*pi/golden_ratio,ticks=[[-pi,pi],None])
both=graphics_array([pccin,wccin])
both.save(frame=true,tick_formatter='latex',fontsize=12,filename='junk_sage.pdf')
