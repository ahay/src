#!/usr/bin/gnuplot
set terminal postscript eps size 3,3 enhanced font "Times-Italic"
set output 
#set output 'opwd.eps'
#unset border
set size square
unset key
set grid 
set ytics ("x_1-{/Symbol D}x_1" -1,"x_1" 0,"x_1+{/Symbol D}x_1" 1)
set xtics ("x_2-{/Symbol D}x_2" -1,"x_2" 0,"x_2+{/Symbol D}x_2" 1)
set parametric
theta=pi/6
set xrange [-2:2]
set yrange [-2:2]
set trange [-1:1]
set label 1 "O" at -0.11,0.08
set label 3 "A" at 1.1,tan(theta)
set label 4 "B" at cos(theta)-0.2,sin(theta)
set label 5 "{/Symbol \161}" at 0.25,0.08
set label 11 at 0,0 point pt 7
set label 12 at 0,1 point pt 7 
set label 13 at 0,-1 point pt 7 
set label 14 at -1,0 point pt 7 
set label 15 at -1,1 point pt 7 
set label 16 at -1,-1 point pt 7 
set label 17 at 1,0 point pt 7 
set label 18 at 1,1 point pt 7 
set label 19 at 1,-1 point pt 7 
set label 41 at 1,tan(theta) point pt 7 ps 0.7
set label 42 at cos(theta),sin(theta) point pt 7 ps 0.7
func(x)=tan(theta)*x
dt=0.01
f0=pi*30*dt
ricker(x)=(1-2*f0*f0*x*x)*exp(-f0*f0*x*x)
x1(x)=cos(x*2*pi)
y1(x)=sin(x*2*pi)
x2(x)=x*1.7
y2(x)=func(x2(x))
x3(x)=0.2*cos(0.5*theta*(x+1))
y3(x)=0.2*sin(0.5*theta*(x+1))
y4(x)=2*x
x4(x)=0.1*ricker(y4(x))
x5(x)=0.1*ricker(y4(x)-func(1))+1
x6(x)=0.1*ricker(y4(x)+func(1))-1
plot x1(t),y1(t), x2(t),y2(t) ls 1, x3(t),y3(t) ls 1, x4(t),y4(t) ls 1, x5(t),y4(t) ls 1, x6(t),y4(t) ls 1

