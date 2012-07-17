import m8r

# Download data
m8r.Fetch('imgs','mona.img')
fmona=open('mona.txt','wt')
fmona.write('n1=512 n2=513 in=mona.img data_format=native_uchar\n')
fmona.close()

mona = m8r.dd(type='float').transp['mona.txt']
mona.grey(title='Mona Lisa',allpos=1,color='b',screenratio=1,wantaxis=0).show()

# Edge preserving smoothing 
mona2 = m8r.impl2(rect1=80,rect2=80,tau=1)[mona]
mona2.grey(title='Smoothed',allpos=1,color='b',screenratio=1,wantaxis=0).show()

