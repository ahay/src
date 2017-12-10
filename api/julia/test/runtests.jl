#!/usr/bin/env julia
using Base.Test

Base.append!(ARGS, ["int1=1", "float1=1e-99", "str1=ḉ", "bool1=n"])

using m8r

println("RSFROOT")
@test m8r.RSFROOT ≠ nothing

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

println("input int")
run(pipeline(pipeline(`sfspike n1=2 k1=1,2,2
                               n2=3 k2=1,2,3
                               nsp=3 mag=1,4,2`,
                      `sfdd type=int`),
             stdout="test_inp_int.rsf"))
inp = m8r.input("test_inp_int.rsf")
@test m8r.size(inp) == (2, 3)
@test m8r.gettype(inp) == 3
dat, = m8r.read(inp)
@test dat == Int32[1 0 0; 0 4 2]
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
run(`sfrm test_inp_int.rsf`)

println("output int")
out = m8r.output("test_out_int.rsf")
m8r.setformat(out, "int")
@test m8r.gettype(out) == 3
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

m8r.intwrite(Int32[1; 2], Int32[m8r.leftsize(out, 0)][], out)

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
m8r.close(out)

stdout = STDOUT
(rout, wout) = redirect_stdout()
run(pipeline(`sfdisfil`, stdin="test_out_int.rsf", stdout=wout))
data = convert(String, readavailable(rout))
close(rout)
redirect_stdout(stdout)
@test data == "   0:    1    2 \n"
run(`sfrm test_out_int.rsf`)


println("input float")
run(pipeline(`sfspike n1=2 k1=1,2,2
                      n2=3 k2=1,2,3
                      nsp=3 mag=1,4,2 out=stdout`,
             stdout="test_inp_float.rsf"))
inp = m8r.input("test_inp_float.rsf")
@test m8r.size(inp) == (2, 3)
@test m8r.gettype(inp) == 4
dat, = m8r.read(inp)
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
run(`sfrm test_inp_float.rsf`)

println("output float")
out = m8r.output("test_out_float.rsf")
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

m8r.floatwrite(Float32[1.5; 2.5], Int32[m8r.leftsize(out, 0)][], out)

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
m8r.close(out)

stdout = STDOUT
(rout, wout) = redirect_stdout()
run(pipeline(`sfdisfil`, stdin="test_out_float.rsf", stdout=wout))
data = convert(String, readavailable(rout))
close(rout)
redirect_stdout(stdout)
@test data == "   0:           1.5          2.5\n"
run(`sfrm test_out_float.rsf`)

println("input complex")
open("test_inp_complex.rsf", "w") do rsf_f
    run(pipeline(`sfspike n1=2 k1=1,2,2
                               n2=3 k2=1,2,3
                               nsp=3 mag=1,4,2`,
                       `sfrtoc`,
                       `sfmath output='input + I' out=stdout`, rsf_f))
end
inp = m8r.input("test_inp_complex.rsf")
@test m8r.size(inp) == (2, 3)
@test m8r.gettype(inp) == 5
dat, = m8r.read(inp)
@test dat == Complex64[1+im im im; im 4+im 2+im]

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
run(`sfrm test_inp_complex.rsf`)

println("output complex")
out = m8r.output("test_out_complex.rsf")
m8r.setformat(out, "complex")
@test m8r.gettype(out) == 5
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

m8r.complexwrite(Complex64[0.5+im; 2+im], Int32[m8r.leftsize(out, 0)][], out)

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
m8r.close(out)

stdout = STDOUT
(rout, wout) = redirect_stdout()
run(pipeline(`sfdisfil`, stdin="test_out_complex.rsf", stdout=wout))
data = convert(String, readavailable(rout))
close(rout)
redirect_stdout(stdout)
@test data == "   0:        0.5,         1i         2,         1i\n"
run(`sfrm test_out_complex.rsf`)

println("prog")
println("    read")
dat, n, d, o, l, u = sfspike(n1=4, n2=2, nsp=2, k1=(1,2), mag=(1,3)) |>
                     x -> sfwindow(x; n1=1, squeeze=false) |>
                     read
@test dat ≈ [1. 1.]
@test n == [1, 2]
@test d ≈ [0.004, 0.1]
@test o ≈ [0, 0]
@test l == ["Time", "Distance"]
@test u == ["s", "km"]

dat, n, d, o, l, u = sfspike(n1=4, nsp=2, k1=(1,2), mag=(1,3)) |>
                     sfrtoc |>
                     x -> sfmath(x; output="input + I") |>
                     x -> sfwindow(x; n1=3) |>
                     read
@test dat ≈ [1.0+1.0im, 3.0+1.0im, 1.0im]
@test n == [3]
@test d ≈ [0.004]
@test o ≈ [0]
@test l == ["Time"]
@test u == ["s"]

println("    write")
dat, n, d, o, l, u = write([1.1; 0; 0.5], [1 3], [.1 .2 .3], [.4 .5 .5],
                                  ["a" "b" "c"], ["d" "e" "f"]) |>
                     read
@test dat ≈ [1.1 0 0.5]
@test n == [1, 3]
@test d ≈ [0.1, 0.2]
@test o ≈ [0.4, 0.5]
@test l == ["a", "b"]
@test u == ["d", "e"]

dat, n, d, o, l, u = write([im; 0; 0.5], [1 3], [.1 .2 .3], [.4 .5 .5],
                                  ["a" "b" "c"], ["d" "e" "f"]) |>
                     x -> sfadd(x; scale=2) |>
                     read
@test dat ≈ [2im 0 1]
@test n == [1, 3]
@test d ≈ [0.1, 0.2]
@test o ≈ [0.4, 0.5]
@test l == ["a", "b"]
@test u == ["d", "e"]

sfspike(;n1=1) |> x -> write(x, "test_write.rsf")
dat, n, d, o, l, u = read("test_write.rsf")
@test dat ≈ [1.]
@test n == [1]
@test d ≈ [0.004]
@test o ≈ [0]
@test l == ["Time"]
@test u == ["s"]
run(`sfrm test_write.rsf`)

println("all good!")
