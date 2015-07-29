a,vn,vx=var('a,vn,vx')
vbh(a,vn,vx)=1/sqrt(cos(a)^2 + (1/vn - 1/vx)*cos(a)^2*sin(a)^2 + 1/vx*sin(a)^2)
def pbh(vn,vx):
    return parametric_plot([vbh(a,vn,vx)*sin(a),-vbh(a,vn,vx)*cos(a)],(a,-pi/2,pi/2),thickness=2,aspect_ratio=90/162)

vvz(a,vn,vx)=1/(cos(a)*(1 - vn/vx) + sqrt(sin(a)^2*1/vx + cos(a)^2*vn^2/vx^2))
def pvz(vn,vx):
    return parametric_plot([vvz(a,vn,vx)*sin(a),-vvz(a,vn,vx)*cos(a)],(a,-pi/2,pi/2),linestyle='--',color='black')

p = pbh(1.2,1.4)+pvz(1.2,1.4)
p.save(filename='junk_sage.pdf',frame=True,ticks=[[],[]])
