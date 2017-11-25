#!/usr/bin/env julia
using Base.Test

open("test_inp.rsf", "w") do rsf_f
    Base.run(pipeline(`sfspike n1=2 k1=1,2,2
                               n2=3 k2=1,2,3
                               nsp=3 mag=1,4,2 out=stdout`, rsf_f))
end
Base.append!(ARGS, ["int1=1", "float1=1e-99", "str1=ḉ", "bool1=n",
                    "file=test_inp.rsf"])

using m8r
m8r.init()

println("getfloat")
@test m8r.getfloat("float1") ≈ 0
@test m8r.getfloat("float2", 2) ≈  2
@test m8r.getfloat("float3", val=3) ≈ 3
@test m8r.getfloat("float4") ≈ 0

println("getint")
@test m8r.getint("int1") == 1
@test m8r.getint("int2", 2) == 2
@test m8r.getint("int3", val=3) == 3
@test m8r.getint("int4") == 0

println("getstring")
@test m8r.getstring("str1") == "ḉ"
@test m8r.getstring("str2", "2") == "2"
@test m8r.getstring("str3", val="3") == "3"
@test m8r.getstring("str4") == ""

println("getbool")
@test m8r.getbool("bool1") == false
@test m8r.getbool("bool2", false) == false
@test m8r.getbool("bool3", val=false) == false
@test m8r.getbool("bool4") == true


println("input")
inp = m8r.input(m8r.getstring("file"))
@test m8r.shape(inp) == [2, 3]
@test m8r.gettype(inp) == 4
dat = m8r.read(inp)
@test dat == Float32[1 0 0; 0 4 2]

@test m8r.histint(inp, "n1") == 2
@test m8r.histint(inp, "n2") == 3
@test m8r.histfloat(inp, "d1") ≈ 0.004
@test m8r.histfloat(inp, "d2") ≈ 0.1
@test m8r.histfloat(inp, "o1") ≈ 0
@test m8r.histfloat(inp, "o2") ≈ 0
@test m8r.histstring(inp, "label1") == "Time"
@test m8r.histstring(inp, "label2") == "Distance"
@test m8r.histstring(inp, "unit1") == "s"
@test m8r.histstring(inp, "unit2") == "km"

println("output")
out = m8r.output("test_out.rsf")
m8r.putint(out, "n1", 1)
m8r.putint(out, "n2", 2)
m8r.putfloat(out, "d1", 3)
m8r.putfloat(out, "d2", 4)
m8r.putfloat(out, "o1", 5)
m8r.putfloat(out, "o2", 6)
m8r.putstring(out, "label1", "a")
m8r.putstring(out, "label2", "é")
m8r.putstring(out, "unit1", "普通话")
m8r.putstring(out, "unit2", "µm")

m8r.floatwrite(Float32[1; 2], Int32[m8r.leftsize(out, 0)][1], out)

@test m8r.histint(out, "n1") == 1
@test m8r.histint(out, "n2") == 2
@test m8r.histfloat(out, "d1") ≈ 3
@test m8r.histfloat(out, "d2") ≈ 4
@test m8r.histfloat(out, "o1") ≈ 5
@test m8r.histfloat(out, "o2") ≈ 6
@test m8r.histstring(out, "label1") == "a"
@test m8r.histstring(out, "label2") == "é"
@test m8r.histstring(out, "unit1") == "普通话"
@test m8r.histstring(out, "unit2") == "µm"

Base.run(`sfrm test_inp.rsf`)
Base.run(`sfrm test_out.rsf`)
println("all good!")
