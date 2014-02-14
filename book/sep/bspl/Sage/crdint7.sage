z1(x)=(2416 - 1680*x^2 + 560*x^4 - 140*x^6 + 35*x^7)/5040
z2(x)=(2472 - 392*x - 504*x^2 - 1960*x^3 + 2520*x^4 - 1176*x^5 + 252*x^6 - 21*x^7)/5040
z3(x)=(-1112 + 12152*x - 19320*x^2 + 13720*x^3 - 5320*x^4 + 1176*x^5 - 140*x^6 + 7*x^7)/5040
z4(x)=(4-x)^7/5040
spl(x)=(heaviside(x+1)-heaviside(x-1))*z1(abs(x))+(heaviside(x-1)-heaviside(x-2))*z2(x)+(heaviside(x-2)-heaviside(x-3))*z3(x)+(heaviside(x-3)-heaviside(x-4))*z4(x)+(heaviside(x+2)-heaviside(x+1))*z2(-x)+(heaviside(x+3)-heaviside(x+2))*z3(-x)+(heaviside(x+4)-heaviside(x+3))*z4(-x)
a=0.535281
b=0.122555
c=0.00914759
d=0.330597
k=var('k')
f(z)=sum((-1)^k*(a*z)^k,k,0,100)*sum((-1)^k*(a/z)^k,k,0,100)*sum((-1)^k*(b*z)^k,k,0,100)*sum((-1)^k*(b/z)^k,k,0,100)*sum((-1)^k*(c*z)^k,k,0,100)*sum((-1)^k*(c/z)^k,k,0,100)/d
t=taylor(f(x),x,0,100)
g(x)=sum([t.coefficient(x^k)*spl(x-k) for k in range(-10,0)+range(1,11)])+spl(x)*4.9647378749
pspl=plot(g(x),x,-6,6,thickness=3,aspect_ratio=10/golden_ratio)
splf(w)=(2/w*sin(w/2))^8/(z1(0) + 2*z2(1)*cos(w) + 2*z3(2)*cos(2*w) + 2*z4(3)*cos(3*w))
wspl=plot(splf(w),w,-2*pi,2*pi,thickness=3,color='red',aspect_ratio=4*pi/golden_ratio,ticks=[[-pi,pi],None])
both=graphics_array([pspl,wspl])
both.save(frame=true,tick_formatter='latex',fontsize=12,filename='junk_sage.pdf')