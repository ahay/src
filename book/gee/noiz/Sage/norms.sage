d = [1,1,2,3,5]
m=var('m')

def l2(m,n):
    return sum ([(d[i] - m)^2 for i in range(n)])
def l1(m,n):
    return sum ([abs(d[i] - m) for i in range(n)])
def lp(m,n):
    return sum ([abs(d[i] - m)^(0.1) for i in range(n)])

a = plot(l2(m,1),(m,0,2))+text('$L_2$',(1,1),fontsize=14)
b = plot(l1(m,1),(m,0,2))+text('$L_1$',(1,1),fontsize=14)
c = plot(lp(m,1),(m,0,2),plot_points=200)+text('$L_{0.1}$',(1,1),fontsize=14)

aa = plot(l2(m,5),(m,0,6))
bb = plot(l1(m,5),(m,0,6))
cc = plot(lp(m,5),(m,0,6),plot_points=200)

p = graphics_array([[a,b,c],[aa,bb,cc]])
p.save(filename='junk_sage.pdf')