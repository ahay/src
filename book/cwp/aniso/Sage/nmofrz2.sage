a,vn,vx=var('a,vn,vx')
tbh(x,vx,vn)=sqrt(1 + x^2/vn - x^4/(1 + x^2)*(1/vn - 1/vx))
def pltbh(vx,vn):
    return plot(-tbh(x,vx,vn),(x,-3,3),thickness=2,aspect_ratio=210/162)
		
tvz(x,vx,vn)=(1 - vn/vx) + vn/vx*sqrt(1 + vx*x^2 /vn^2)
def pltvz(vx,vn):
    return plot(-tvz(x,vx,vn),(x,-3,3),linestyle='--',color='black')

p = pltbh(1.2,1.4)+pltvz(1.2,1.4)
p.save(filename='junk_sage.pdf',frame=True,ticks=[[],[]])
