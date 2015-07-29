z1(x)=(2416 - 1680*x^2 + 560*x^4 - 140*x^6 + 35*x^7)/5040
z2(x)=(2472 - 392*x - 504*x^2 - 1960*x^3 + 2520*x^4 - 1176*x^5 + 252*x^6 - 21*x^7)/5040
z3(x)=(-1112 + 12152*x - 19320*x^2 + 13720*x^3 - 5320*x^4 + 1176*x^5 - 140*x^6 + 7*x^7)/5040
z4(x)=(4-x)^7/5040
spl(x)=(heaviside(x+1)-heaviside(x-1))*z1(abs(x))+(heaviside(x-1)-heaviside(x-2))*z2(x)+(heaviside(x-2)-heaviside(x-3))*z3(x)+(heaviside(x-3)-heaviside(x-4))*z4(x)+(heaviside(x+2)-heaviside(x+1))*z2(-x)+(heaviside(x+3)-heaviside(x+2))*z3(-x)+(heaviside(x+4)-heaviside(x+3))*z4(-x)
pspl=plot(spl(x),x,-6,6,thickness=3,aspect_ratio=24/golden_ratio)
splf(w)=(2/w*sin(w/2))^8
wspl=plot(splf(w),w,-2*pi,2*pi,thickness=3,color='red',aspect_ratio=4*pi/golden_ratio,ticks=[[-pi,pi],None])
both=graphics_array([pspl,wspl])
both.save(frame=true,tick_formatter='latex',fontsize=12,filename='junk_sage.pdf')
