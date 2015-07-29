z1(x)=(4-6*x^2+3*x^3)/6
z2(x)=(2-x)^3/6
spl(x)=(heaviside(x+1)-heaviside(x-1))*z1(abs(x))+(heaviside(x-1)-heaviside(x-2))*z2(x)+(heaviside(x+2)-heaviside(x+1))*z2(-x)
w=var('w')
splf(w)=2*(integrate(cos(w*x)*z1(x),x,0,1)+integrate(cos(w*x)*z2(x),x,1,2))/(z1(0) + 2*z1(1)*cos(w))
a = n(2 - sqrt(3))
c = n((2+ sqrt(3))/6)
k,m,z=var('k,m,z')
f(z)=sum((-1)^k*(a*z)^k,k,0,100)*sum((-1)^k*(a/z)^k,k,0,100)/c
t=taylor(f(x),x,0,100)
g(x)=sum([t.coefficient(x^k)*spl(x-k) for k in range(-10,0)+range(1,11)])+spl(x)*1.73205080757
pspl=plot(g(x),x,-6,6,thickness=3,aspect_ratio=10/golden_ratio)
wspl=plot(splf(w),w,-2*pi,2*pi,thickness=3,color='red',aspect_ratio=4*pi/golden_ratio,ticks=[[-pi,pi],None])
both=graphics_array([pspl,wspl])
both.save(frame=true,tick_formatter='latex',fontsize=12,filename='junk_sage.pdf')