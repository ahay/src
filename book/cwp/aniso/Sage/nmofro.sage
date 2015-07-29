a,v,vn,vx=var('a,v,vn,vx')
vtt(a,vn,vx)=1/sqrt(cos(a)^2 + (1/vn - 1/vx)*cos(a)^2*sin(a)^2 + 1/vx*sin(a)^2)
def ptt(vn,vx,c):
    return parametric_plot([vtt(a,vn,vx)*sin(a),-vtt(a,vn,vx)*cos(a)],(a,-pi/2,pi/2),thickness=3)+line([(-vx-0.01,0),(vx+0.01,0)],color='black')+line([(0,0),(0,-vx-0.01)],color='black')+text(c,(0,0.1),fontsize=14)
def izo(v):
    return parametric_plot([sqrt(v)*sin(a),-sqrt(v)*cos(a)],(a,-pi/2,pi/2),linestyle='--',color='black')
p=graphics_array([[ptt(1,1,'a'),ptt(1.4,1.4,'b')+izo(1)+izo(1.4)],[ptt(1,1.4,'c')+izo(1)+izo(1.4),ptt(0.6,1.4,'d')+izo(1)+izo(1.4)]])
p.save(filename='junk_sage.pdf',axes=False)