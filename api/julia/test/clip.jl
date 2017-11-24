#!/usr/bin/env julia

using m8r

m8r.init()
inp = m8r.input("in")
out = m8r.output("out")

n1 = m8r.histint(inp,"n1")
n2 = m8r.leftsize(inp,1)

clip = m8r.getfloat("clip")

trace = Array{Float32}(n1)

for i2 in 1:n2
    m8r.floatread(trace,n1,inp)
    trace = clamp(trace,-clip,clip)
    m8r.floatwrite(trace,n1,out)
end



