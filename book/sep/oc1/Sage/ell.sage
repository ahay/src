hh(x)=-sqrt(1-x^2/2)
ph=plot(hh(x),(x,-0.999999*sqrt(2),0.999999*sqrt(2)),thickness=3,color='green')
h,t=var('h,t')
r(h)=(0.01+1-h+sign(h-1)*sqrt((0.01-1-h)^2-4*h))/0.1
def pr(s):
    return parametric_plot((0.1+sqrt(s)+(r(s)-0.1-sqrt(s))*t,t*hh(r(s))),(t,0,1))+parametric_plot((0.1-sqrt(s)+(r(s)-0.1+sqrt(s))*t,t*hh(r(s))),(t,0,1))
p=pr(0.16)+pr(0.64)+pr(1.44)+ph
p.save(ticks=[None,[]],filename='junk_sage.pdf')
