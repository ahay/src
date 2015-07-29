sinc(x) = sin(pi*x)/(pi*x)
n=var('n')
wind(x,n) = (pi*x)/(2*n*tan((pi*x)/(2*n)))
muir(x,n) = sin(pi*x)/(2*n*tan((pi*x)/(2*n)))

def ps(n):
    return plot(sinc(x),(x,-n,n),thickness=2)
def pw(n):
    return plot(wind(x,n),(x,-n,n),thickness=2,color='green')
def pm(n):
    return plot(muir(x,n),(x,-n,n),thickness=2,color='red')

p=graphics_array([[ps(2),pw(2),pm(2)],
		  [ps(6),pw(6),pm(6)]])
p.save(filename='junk_sage.pdf',ymin=-0.25,ymax=1.05,frame=True)