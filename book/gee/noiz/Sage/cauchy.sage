d=[1,1,2,3,5]
m=var('m')

def A(r,n,m):
    return sum ([log(1 + (d[i] - m)^2/r^2) for i in range(n)])

def a1(r):
    return plot(A(r,1,m),(m,0,2),ymin=0)
def an(r):
    return plot(A(r,5,m),(m,0,6),ymin=0)

p = graphics_array([[a1(2),a1(1),a1(0.2)],[an(2),an(1),an(0.2)]])
p.save(filename='junk_sage.pdf')