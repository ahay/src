# ------------------------------------------------------------
# strings

a='Golden'
len(a)
a[0]
a[0:4]

b = a + ' ' + 'workshop'
print b

c = b + 2018
c = b + ' '+str(2018)
print c

# ------------------------------------------------------------
# lists

d = ['Golden', 'workshop']
len(d)
print d[0]
print d[1]

d.append('2018')
print d

# ------------------------------------------------------------
# tuple = a sequence of immutable Python objects. 
t = ('Golden', 'workshop')

t = t + (2018,)
print t

# ------------------------------------------------------------
# dictionaries

e = {'what':'workshop','where':'Golden','when':2018}
print e
print e['where']+' '+e['what']+' '+str(e['when'])

f = dict(what='workshop',where='Melbourne',when=2013)
print f
print f['where']+' '+f['what']+' '+str(f['when'])

# ------------------------------------------------------------
# loops

for i in range(len(a)):
    print a[i]

for i in range(len(d)):
    print d[i]

for i in t:
    print i
    
for key in e.keys():
    print key,e[key]


# ------------------------------------------------------------
# conditional statements

for k in range(5):
    if k < 2:
        print k,'<2'
    else:
        print k,'>=2'

try:
    b+2014
except:
    print "error!"

# ------------------------------------------------------------
# functions

def m8rschool(year):
    workshops = dict(Golden=2018,Melbourne=2013)
    for key in workshops.keys():
        if workshops[key]==year:
            return key            

print m8rschool(2018)


def increment(a,b=5):
    return a+b

# ------------------------------------------------------------
# modules
import math
x = math.sqrt(increment(4))
print x
