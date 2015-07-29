s,r=var('s,r')
f(r)=sqrt((r-0.5)*(r+0.5))
bf(r)=sqrt((0.5-r)*(r+0.5))
t(s,r)=sqrt((r-s)^2-(f(s)+f(r))^2)
bt(s,r)=sqrt((r-s)^2+(bf(s)+bf(r))^2)
def pt(s):
    return plot(abs(t(s,r)),(r,-0.5,0.5))
p1=add(map(pt,[-0.4,0.2,0.0,0.2,0.4]))
def ptm(s):
    return plot(abs(bt(s,r)),(r,-1.5,-0.5))
p2=add(map(ptm,[0.6,0.8,1.0]))
def ptp(s):
    return plot(abs(bt(s,r)),(r,0.5,1.5))
p3=add(map(ptp,[-0.6,-0.8,-1.0]))
p=p1+p2+p3+text("-0.4",(-0.625,0.3),fontsize=14)+text("-0.2",(-0.625,0.55),fontsize=14)+text("0",(-0.575,0.7),fontsize=14)+text("0.2",(-0.625,0.85),fontsize=14)+text("0.4",(-0.625,0.95),fontsize=14)+text("0.6",(-0.4,1.075),fontsize=14)+text("0.8",(-0.4,1.15),fontsize=14)+text("1.0",(-0.4,1.225),fontsize=14)+text("-0.6",(0.375,1.075),fontsize=14)+text("-0.8",(0.375,1.15),fontsize=14)+text("-1.0",(0.375,1.225),fontsize=14)
p.axes_range(ymin=1.5,ymax=0)
p.save(axes_labels=['$r$',''],fontsize=14,figsize=[12,6],filename='junk_sage.pdf')