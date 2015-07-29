# ------------------------------------------------------------
# strings

a='StPetersburg'
len(a)
a[0]
a[4:7]

b=a+' '+'workshop'
print b

c=b+2014
c=b+' '+str(2014)
print c

# ------------------------------------------------------------
# lists

d = ['StPetersburg', 'workshop']
len(d)
print d[0]
print d[1]

d.append('2014')
print d

# ------------------------------------------------------------
# tuple = a sequence of immutable Python objects. 
t = ('StPetersburg', 'workshop')

t = t + (2014,)
print t

# ------------------------------------------------------------
# dictionaries

e={'what':'workshop','where':'StPetersburg','when':2014}
print e
print e['where']+' '+e['what']+' '+str(e['when'])

f=dict(what='workshop',where='Melbourne',when=2013)
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
    workshops=dict(StPetersburg=2014,Melbourne=2013)
    for key in workshops.keys():
        if workshops[key]==year:
            return key            

print m8rschool(2014)


def increment(a,b=5):
    return a+b

# ------------------------------------------------------------
# modules
import math
x=math.sqrt(increment(4))
print x
