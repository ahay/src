#!/usr/bin/env julia

using RegularlySampledFormat

inp = Input("in")
out = Output("out")

n1 = getint(inp, "n1")
n2 = leftsize(inp, 1)

clip = getfloat("clip")
if clip == :none
   error("Need clip=")
end

trace = Array{Float32}(undef, n1)

for i2 in 1:n2
    dataread!(inp, trace)
    clamp!(trace, -clip, clip)
    datawrite(out, trace)
end

fileclose(inp)
fileclose(out)

