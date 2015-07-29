x,h,t=var('x,h,t')
z(x,h,t) = -sqrt((-h + t)*(h + t)*(t - x)*(t + x))/t
a=var('a')
def hyper(a,h,t):
  return plot(-sqrt((-h^2 + t^2)*cos(a)^2 + (x - t*sin(a))^2),(x,-2*t,2*t),color='black')
hp = sum([hyper(a*pi/20,0.8,1.2) for a in range(-20,21)])
zp = plot(z(x,0.8,1.2), (x, -1.2, 1.2),thickness=4,color='red') 
yp = plot(-sqrt(1.2^2 - 0.8^2)*sqrt((0.8 - x)*(0.8 + x))/0.8,(x, -0.8, 0.8),thickness=4,color='green') 
ap = arrow((-1.5,0),(1.5,0))+line([(0,0),(0,-2.1)])
tp = text('$-h$',(-0.8,0.05),fontsize=16)+text('$h$',(0.8,0.05),fontsize=16)+text('$y$',(1.55,0),fontsize=16)
ep=ellipse((-0.8,0),0.02,0.01,fill=True)+ellipse((0.8,0),0.02,0.01,fill=True)
p=hp+zp+yp+ap+tp+ep
p.save(filename='junk_sage.pdf',xmin=-1.5,xmax=1.6,ymin=-1,ymax=0.2,aspect_ratio='automatic',axes=False)
