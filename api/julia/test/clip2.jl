#!/usr/bin/env julia

using RegularlySampledFormat

inp = Input("in")
out = Output("out")

n1 = getint(inp, "n1")
n2 = leftsize(inp, 1)

clip = getfloat("clip")

trace = Array{Float32}(undef, n1)

for i2 in 1:n2
    floatread(inp, trace)
    clamp!(trace, -clip, clip)
    floatwrite(out, trace)
end

fileclose(inp)
fileclose(out)

