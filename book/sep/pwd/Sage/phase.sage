p,z,w=var('p,z,w')
m(p)=1/12*(1+p)*(2+p)
m1(p)=1/1680*(1+p)*(2+p)*(3+p)*(4+p)
m2(p)=1/420*(4-p)*(2+p)*(3+p)*(4+p)
f3(p,z)=m(p)*z + (1-m(p)-m(-p)) + m(-p)/z
f5(p,z)=m1(p)*z^2 + m2(p)*z + (1-m1(p)-m1(-p)-m2(p)-m2(-p)) + m2(-p)/z + m1(-p)/z^2
a3(p,z)=f3(p,z)/f3(p,1/z)
a5(p,z)=f5(p,z)/f5(p,1/z)
p1=plot(w*1/2,(w,0,pi),legend_label='exact',color='black')+plot(arg(a3(1/2,exp(I*w))),(w,0,pi),legend_label='3-point',linestyle='--',color='green')+plot(arg(a5(1/2,exp(I*w))),(w,0,pi),legend_label='5-point',linestyle=':')
p1.set_legend_options(loc='upper center')
p2=plot(w*4/5,(w,0,pi),color='black')+plot(arg(a3(4/5,exp(I*w))),(w,0,pi),linestyle='--',color='green')+plot(arg(a5(4/5,exp(I*w))),(w,0,pi),linestyle=':')
p=graphics_array([p1,p2])
p.save(frame=True,axes_labels=['frequency','phase'],figsize=[12,6],fontsize=14,filename='junk_sage.pdf')