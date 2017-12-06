__precompile__()
"""
    m8r.jl

Julia interface to Madagascar
"""
module m8r

import Base.read
import Base.write

export read,
       write

if haskey(ENV, "RSFROOT")
    RSFROOT = ENV["RSFROOT"]
else
    RSFROOT = nothing
end

immutable File
    tag::String
    rsf::Ptr{UInt8}
end

function __init__()
    src = Base.source_path()
    argv = src == nothing ? ["julia"] : [basename(src)]
    if isempty(ARGS)
        push!(argv, "-")
    else
        append!(argv,ARGS)
    end
    ccall((:sf_init,"libdrsf"),Void,(Int32,Ptr{Ptr{UInt8}}),length(argv),argv)
end

function input(tag::String)
    if tag ≠ "in"
        if !isfile(tag)
            throw("SystemError: unable to read file $tag")
        end
    end
    rsf = ccall((:sf_input,"libdrsf"),Ptr{UInt8},(Ptr{UInt8},),tag)
    File(tag,rsf)
end

function output(tag::String)
    rsf = ccall((:sf_output,"libdrsf"),Ptr{UInt8},(Ptr{UInt8},),tag)
    File(tag,rsf)
end

function setformat(file::File,format::String)
    ccall((:sf_setformat,"libdrsf"),Void,(Ptr{UInt8},Ptr{UInt8}),file.rsf,format)
end

function histint(file::File,name::String)
    val = Cint[0]
    ccall((:sf_histint,"libdrsf"),Bool,(Ptr{UInt8},Ptr{UInt8},Ptr{Cint}),file.rsf,name,val)
    return val[]
end

function histfloat(file::File,name::String)
    val = Cfloat[0]
    ccall((:sf_histfloat,"libdrsf"),Bool,(Ptr{UInt8},Ptr{UInt8},Ptr{Cfloat}),file.rsf,name,val)
    return val[]
end

function histstring(file::File,name::String)
    val = ccall((:sf_histstring,"libdrsf"),Ptr{Cchar},(Ptr{UInt8},Ptr{UInt8}),file.rsf,name)
    if val == C_NULL
        return ""
    end
    return unsafe_string(val)
end

function getint(name::String, val::Int)
    val = Cint[val]
    ccall((:sf_getint,"libdrsf"),Bool,(Ptr{UInt8},Ptr{Cint}),name,val)
    return val[]
end
getint(name; val::Int = 0) = getint(name, val)

function getfloat(name::String, val::Real)
    val = Cfloat[val]
    ccall((:sf_getfloat,"libdrsf"),Bool,(Ptr{UInt8},Ptr{Cfloat}),name,val)
    return val[]
end
getfloat(name::String; val::Real = 0) = getfloat(name, val)

function getstring(name::String, val::String)
    v = ccall((:sf_getstring,"libdrsf"),Ptr{Cchar},(Ptr{UInt8},),name)
    if v == C_NULL
        return val
    end
    return unsafe_string(v)
end
getstring(name::String; val::String = "") = getstring(name, val)

function getbool(name::String, val::Bool)
    val = Bool[val]
    ccall((:sf_getbool,"libdrsf"),Bool,(Ptr{UInt8},Ptr{Bool}),name,val)
    return val[]
end
getbool(name::String; val::Bool = true) = getbool(name, val)

function gettype(file::File)
    return ccall((:sf_gettype,"libdrsf"),Cuint,(Ptr{UInt8},),file.rsf) + 1
end

function leftsize(file::File,dim::Int)
    ccall((:sf_leftsize,"libdrsf"),Culonglong,(Ptr{UInt8},Cint),file.rsf,dim)
end

function floatread(arr::Array{Float32,1},size::Int32,file::File)
    ccall((:sf_floatread,"libdrsf"),Void,(Ptr{Cfloat},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function complexread(arr::Array{Complex64,1},size::Int32,file::File)
    ccall((:sf_complexread,"libdrsf"),Void,(Ptr{Complex64},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function floatwrite(arr::Array{Float32,1},size::Int32,file::File)
    ccall((:sf_floatwrite,"libdrsf"),Void,(Ptr{Cfloat},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function complexwrite(arr::Array{Complex64,1},size::Int32,file::File)
    ccall((:sf_complexwrite,"libdrsf"),Void,(Ptr{Complex64},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function putint(file::File,name::String,val::Int)
    val = convert(Cint, val)
    ccall((:sf_putint,"libdrsf"),Void,(Ptr{UInt8},Ptr{UInt8},Cint),file.rsf,name,val)
end

function putfloat(file::File,name::String,val::Real)
    val = convert(Cfloat, val)
    ccall((:sf_putfloat,"libdrsf"),Void,(Ptr{UInt8},Ptr{UInt8},Cfloat),file.rsf,name,val)
end

function putstring(file::File,name::String,val::String)
    ccall((:sf_putstring,"libdrsf"),Void,(Ptr{UInt8},Ptr{UInt8},Ptr{UInt8}),file.rsf,name,val)
end

function close(file::File)
    ccall((:sf_fileclose,"libdrsf"), Void, (Ptr{UInt8},), file.rsf)
end

"""
    m8r.size(file::m8r.File) -> Tuple

The size of `file`, an Int32 array representing the length of each of its
dimensions.

# Examples

```julia-repl
julia> open("spike.rsf", "w") do rsf_f
run(pipeline(`sfspike n1=2 n2=3`, rsf_f))
end

julia> inp = m8r.input("spike.rsf")

julia> m8r.size(inp)
(2, 3)
```
"""
function size(file::File)
    size = leftsize(file, 0)
    dim = 1
    n = histint(file, string("n", dim))
    s = Int32[n]
    size /= n
    dim += 1
    while size > 1
        n = histint(file, string("n", dim))
        push!(s, n)
        size /= n
        dim += 1
    end
    return Tuple(s)
end

"""
    read(file)

Reads RSF `file`, returning its contents and header. `file` may be file handle 
(`m8r.File`), filename (`m8r.File.tag`), or `NTuple{2, Base.PipeEndpoint}`. This
last option is useful for reading pipes from `sf` commands (see last example).

# Examples

## Reading file handle
```julia-repl
julia> open("spike.rsf", "w") do rsf_f
run(pipeline(`sfspike n1=2 n2=3`, rsf_f))
end

julia> inp = input("spike.rsf")

julia> dat, n, d, o, l, u = read(inp)
(Float32[1.0 1.0 1.0; 1.0 1.0 1.0], [2, 3], Float32[0.004, 0.1], Float32[0.0, 0.0], String["Time", "Distance"], String["s", "km"])
```

## Reading file name
```julia-repl
julia> read("spike.rsf")
(Float32[1.0 1.0 1.0; 1.0 1.0 1.0], [2, 3], Float32[0.004, 0.1], Float32[0.0, 0.0], String["Time", "Distance"], String["s", "km"])
```

## Reading pipe
```julia-repl
julia> sfspike(;n1=1, n2=2) |> read
(Float32[1.0 1.0 1.0; 1.0 1.0 1.0], [2, 3], Float32[0.004, 0.1], Float32[0.0, 0.0], String["Time", "Distance"], String["s", "km"])
```
"""
function read(file::File)
    t = [
         UInt8, # SF_UCHAR
         UInt8, # SF_CHAR
         Int, # SF_INT
         Float32, # SF_FLOAT
         Complex64, # SF_COMPLEX
         Int16, # SF_SHORT
         Float64, # SF_DOUBLE
         Int, # SF_LONG
        ]
    n = Int[i for i in size(file)]
    sz::Int32 = prod(n)
    t_idx = gettype(file)
    data = zeros(t[t_idx], sz)
    if t_idx == 4
        floatread(data, sz, file)
    elseif t_idx == 5
        complexread(data, sz, file)
    else
        throw("Can only read Float32 and Complex64 (not implemented)")
    end

    data = reshape(data, n...)
    n = Int[i for i in size(file)]
    d = Float32[]
    o = Float32[]
    l = String[]
    u = String[]
    for i in 1:length(n)
        append!(d, [histfloat(file, "d"*dec(i))])
        append!(o, [histfloat(file, "o"*dec(i))])
        append!(l, [histstring(file, "label"*dec(i))])
        append!(u, [histstring(file, "unit"*dec(i))])
    end
    return data, n, d, o, l, u
end

read(name::String) = read(input(name))

function read(stdin::NTuple{2, Base.PipeEndpoint})
    rin, win = stdin
    flush(win)
    old_stdin = STDIN
    redirect_stdin(rin)
    data = read("in")
    redirect_stdin(old_stdin)
    return data
end

"""
    write(file, dat, n, d, o, label, unit)

Write RSF `file`. `file` may be file handle (`m8r.File`), filename
(`m8r.File.tag`) or absent:

    write(dat, n, d, o, label, unit) -> NTuple{2, Base.PipeEndpoint}

This last option is useful for writing pipes to `sf` commands (see last
example). However, it must be noted that in those cases `write` can be omitted.

In all methods, `n`, `d`, `o`, `label`, and `unit` are optional. If given, they
should be of type `AbstractArray`.

Finally, one may write from a pipe to a file with:

    write(stdin::NTuple{2, Base.PipeEndpoint}, file::String)

!!! warning "Writing to file handles"

    Do *not* use supply the file as an `m8r.File` type unless you know exactly
    what you are doing. Because of how Madagascar is set up, calling m8r.output
    sets the filetype to whatever was used for the last m8r.input("in") call. If
    this function has not yet been called, it will default to `Float32`.
    Therefore, it is *impossible* to create a `Complex` file afterwards. Other
    `write` methods overcome this limitation by writing dummy files to dummy
    pipes to "trick" Madagascar into switching file types.

    In addition, do not try to write to a filehandle which has alread been
    written to with `write`, as `write` closes the file. Doing so will cause a
    segfault.

# Examples

## Write file by name
```julia-repl
julia> write("spike.rsf", [1., 2.])

julia> read("spike.rsf")
(Float32[1.0, 2.0], [2], Float32[1.0], Float32[0.0], String[""], String[""])
```

## Write to pipe
```julia-repl
julia> write([1. im]) |> sfreal |> read
(Complex{Float32}[1.0+0.0im 0.0+1.0im], [1, 2], Float32[1.0, 1.0], Float32[0.0, 0.0], String["", ""], String["", ""])
```

## Write from pipe
```julia-repl
julia> sfspike(;n1=1) |> x -> write(x, "spike.rsf")

julia> read("spike.rsf")
(Float32[1.0, 2.0], [2], Float32[1.0], Float32[0.0], String[""], String[""])
```

## Writing file handle (avoid this!)
```julia-repl
julia> out = output("spike.rsf")

julia> write(out, [1., 2]) # write(out, [im, 2]) will not work due to warning

julia> read(out.tag)
(Float32[1.0, 2.0], [2], Float32[1.0], Float32[0.0], String[""], String[""])
```
"""
function write(file::File, dat::AbstractArray, n=nothing, d=nothing,
                  o=nothing, l=nothing, u=nothing)
    if n == nothing
        n = size(dat)
    end
    dim = length(n)
    d = d == nothing ? [1 for i in 1:dim] : d
    o = o == nothing ? [0 for i in 1:dim] : o
    l = l == nothing ? ["" for i in 1:dim] : l
    u = u == nothing ? ["" for i in 1:dim] : u
    for i in 1:dim
        putint(file, "n$i", n[i])
        putfloat(file, "d$i", d[i])
        putfloat(file, "o$i", o[i])
        putstring(file, "label$i", l[i])
        putstring(file, "unit$i", u[i])
    end

    if eltype(dat) <: AbstractFloat
        floatwrite(Array{Float32}(vec(dat)), Int32[leftsize(file, 0)][], file)
    elseif eltype(dat) <: Complex
        complexwrite(Array{Complex64}(vec(dat)), Int32[leftsize(file, 0)][],
                     file)
    else
        throw("Can only write Float32 and Complex64 (not implemented)")
    end
    close(file)
end
function write(name::String, dat::AbstractArray, n=nothing, d=nothing,
               o=nothing, l=nothing, u=nothing)
    # This is necessary in case previous command run created a complex file
    spike = joinpath(RSFROOT, "bin", "sfspike")
    old_stdin = STDIN
    (rin, win) = redirect_stdin()
    if  eltype(dat) <: AbstractFloat
        pipe = `$spike n1=1 out=stdout`
    else eltype(dat) <: Complex
        rtoc = joinpath(RSFROOT, "bin", "sfrtoc")
        pipe = pipeline(`$spike n1=1`, `$rtoc out=stdout`)
    end
    p = spawn(pipeline(pipe, stdout=win))
    redirect_stdin(old_stdin)
    Base.wait(p)
    m8r.read((rin, win))

    write(output(name), dat, n, d, o, l, u)
end

function write(dat::AbstractArray, n=nothing, d=nothing, o=nothing,
                  l=nothing, u=nothing)
    if haskey(ENV, "TMPDATAPATH")
        name = joinpath(mktempdir(ENV["TMPDATAPATH"]), "julia.rsf")
    elseif haskey(ENV, "DATAPATH")
        name = joinpath(mktempdir(ENV["DATAPATH"]), "julia.rsf")
    else
        name = joinpath(mktempdir(), "julia.rsf")
    end

    # Slightly roundabout way
    # 1) Write file to disk
    write(name, dat, n, d, o, l, u)

    # 2) Pipe it to dummy sfwindow
    old_stdin = STDIN
    (rin, win) = redirect_stdin()
    progpath = joinpath(RSFROOT, "bin", "sfwindow")
    pipe = `$progpath squeeze=n out=stdout`
    p = spawn(pipeline(pipe, stdin=name, stdout=win))
    redirect_stdin(old_stdin)
    Base.wait(p)

    # 3) Remove temp
    progpath = joinpath(RSFROOT, "bin", "sfrm")
    run(pipeline(`$progpath $name`))
    dir = dirname(name)
    spawn(pipeline(`rmdir $dir`))
    return rin, win
end

function write(stdin::NTuple{2, Base.PipeEndpoint}, name::String)
    dat = read(stdin)
    write(name, dat...)
end

function process_args(;kwargs...)
    args = String[]
    for (key, val) in kwargs
        if typeof(val) <: Tuple
            val = join(["$v" for v in val], ",")
        end
        push!(args, "$key=$val")
    end
    return args
end

if RSFROOT ≠ nothing
    progs = filter(x -> startswith(x, "sf"),
                   readdir(joinpath(RSFROOT, "bin")))
    for (F, S) = [ (Symbol(p), p) for p in progs ]
        @eval export $F
        @eval begin
            progname = $S
"""
    $progname(input) -> NTuple{2, Base.PipeEndpoint}"

Runs RSF program $progname on data provided by `input`. This may be `m8r.File`
or, more commonly `NTuple{2, Base.PipeEndpoint}`. If the program needs no input,
this `input` may be absent.
"""
            function ($F)(stdin::NTuple{2, Base.PipeEndpoint};
                            kwargs...)
                args = process_args(;kwargs...)
                progpath = joinpath(RSFROOT, "bin", $S)
                pipe = `$progpath $args out=stdout`
                rin, win = stdin
                p = spawn(pipeline(pipe, stdin=rin, stdout=win))
                Base.wait(p)
                return rin, win
            end
        end
        @eval function ($F)(;kwargs...)
                args = process_args(;kwargs...)
                progpath = joinpath(RSFROOT, "bin", $S)
                pipe = `$progpath $args out=stdout`
                old_stdin = STDIN
                (rin, win) = redirect_stdin()
                p = spawn(pipeline(pipe, stdout=win))
                redirect_stdin(old_stdin)
                Base.wait(p)
                return rin, win
        end
        @eval function ($F)(file::File; kwargs...)
                args = process_args(;kwargs...)
                progpath = joinpath(RSFROOT, "bin", $S)
                pipe = `$progpath $args out=stdout`
                old_stdin = STDIN
                (rin, win) = redirect_stdin()
                p = spawn(pipeline(pipe, stdin=file.tag, stdout=win))
                redirect_stdin(old_stdin)
                Base.wait(p)
                return rin, win
        end
        @eval function ($F)(dat::AbstractArray, n=nothing, d=nothing, o=nothing,
            l=nothing, u=nothing)
            return write(dat, n, d, o, l, u) |> ($F)
        end
    end
end

end
