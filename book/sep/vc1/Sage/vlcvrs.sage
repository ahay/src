x,v,y=var('x,y,v')
t(x,v)=sqrt(1 + x^2/(1-v^2))
xr(y,v)=y*(1-v^2)
tr(y,v)=sqrt(1 + y^2*(1-v^2))
def vr(y):
    return parametric_plot([xr(y,v),tr(y,v)],(v,0,1.4),linestyle='--',color='red')
vrs=vr(-1)+vr(-0.66)+vr(0.66)+vr(1)
def tt(v):
    return plot(t(x,v),(x,-abs(1-v^2)*1.01,abs(1-v^2)*1.01),thickness=3)
trts=add(map(tt,[0.0,0.6,0.8,1.2,1.4]))
p1=vrs+trts
p1.axes_range(xmin=-1.2,xmax=1.2,ymin=1.6,ymax=0)

pt(x,v)=0.5*x/sqrt(1 - 0.25*v^2)
pxr(y,v)=y*(1 - 0.25*v^2)
ptr(y,v)=0.5*y*sqrt(1 - 0.25*v^2)
def pvr(y):
    return parametric_plot([pxr(y,v),ptr(y,v)],(v,0,2),linestyle='--',color='red')
pvs=add(map(pvr,[0.25,0.50,0.75,1.00]))
def ptt(v):
    return plot(pt(x,v),(x,0,(1-0.25*v^2)*1.1),thickness=3)
prts=add(map(ptt,[0,0.8,1.2,1.6]))
p2=pvs+prts
p2.axes_range(xmin=0,xmax=1.2,ymin=0.6,ymax=0)

p=graphics_array([p1,p2])
p.save(fontsize=14,figsize=[12,6],axes_labels=['x',''],filename='junk_sage.pdf')