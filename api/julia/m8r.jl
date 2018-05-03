__precompile__()
"""
    m8r.jl

Julia interface to Madagascar

Built on the C API, the Julia API provides most lower level functions for
reading, writing and manipulating various types of RSF files. It also provides
higher level read and write functions, documented in `?rsf_read` and
`?rsf_write`.

The Julia API also provides searmless access to Madagascar programs. These can
be accessed as Julia functions. For example, see `?sfwindow`. These functions
can be piped to each other using Julia's currying (function composition)
operator `|>`. They can then be read to memory or written to disk. See examples
below.

# Examples

## Piping to variable
```julia-repl
julia> data = sfspike(n1=2) |> x -> sfwindow(x; n1=1) |> rsf_read

julia> data
(Float32[1.0], [1], Float32[0.004], Float32[0.0], String["Time"], String["s"])
```

## Piping to disk
```julia-repl
julia> sfspike(n1=2) |> x -> sfwindow(x; n1=1) |> x -> rsf_write("spike.rsf", x)
```
Alternatively, and this is especially useful when the contents of the output are too large to be held in RAM, the line above can be replaced with
```julia-repl
julia> sfspike(n1=2) |> x -> sfwindow(x; to_file="spike.rsf", n1=1)
```
"""
module m8r

export rsf_read,
       rsf_write

if haskey(ENV, "RSFROOT")
    RSFROOT = ENV["RSFROOT"]
else
    RSFROOT = nothing
end

#if Libdl.find_library("libdrsf") == ""
    #push!(Libdl.DL_LOAD_PATH, joinpath(RSFROOT, "lib"))
    #if Libdl.find_library("libdrsf") == ""
        #throw("Cannot find C api. Make sure RSFROOT/lib is in LD_LIBRARY_PATH")
    #end
#end

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

function gettype(file::File)
    return ccall((:sf_gettype,"libdrsf"),Cuint,(Ptr{UInt8},),file.rsf) + 1
end

function getform(file::File)
    return ccall((:sf_getform,"libdrsf"),Cuint,(Ptr{UInt8},),file.rsf) + 1
end

function esize(file::File)
    return ccall((:sf_esize,"libdrsf"),Csize_t,(Ptr{UInt8},),file.rsf)
end

function setformat(file::File,format::String)
    ccall((:sf_setformat,"libdrsf"),Void,(Ptr{UInt8},Ptr{UInt8}),file.rsf,format)
end

function histint(file::File,name::String)
    val = Cint[0]
    ccall((:sf_histint,"libdrsf"),Bool,(Ptr{UInt8},Ptr{UInt8},Ref{Cint}),file.rsf,name,val)
    return convert(Int, val[])
end

function histfloat(file::File,name::String)
    val = Cfloat[0]
    ccall((:sf_histfloat,"libdrsf"),Bool,(Ptr{UInt8},Ptr{UInt8},Ref{Cfloat}),file.rsf,name,val)
    return convert(Float32, val[])
end

function histstring(file::File,name::String)
    val = ccall((:sf_histstring,"libdrsf"),Ptr{Cchar},(Ptr{UInt8},Ptr{UInt8}),file.rsf,name)
    if val == C_NULL
        return ""
    end
    return unsafe_string(val)
end

function getint(name::String, val::Integer)
    val = Cint[val]
    ccall((:sf_getint,"libdrsf"),Bool,(Ptr{UInt8},Ref{Cint}),name,val)
    return convert(Int, val[])
end
getint(name::String; val::Integer = 0) = getint(name, val)

function getfloat(name::String, val::Real)
    val = Cfloat[val]
    ccall((:sf_getfloat,"libdrsf"),Bool,(Ptr{UInt8},Ref{Cfloat}),name,val)
    return convert(Float32, val[])
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
    ccall((:sf_getbool,"libdrsf"),Bool,(Ptr{UInt8},Ref{Bool}),name,val)
    return val[]
end
getbool(name::String; val::Bool = true) = getbool(name, val)

function leftsize(file::File,dim::Integer)
    dim::Cint = dim
    ccall((:sf_leftsize,"libdrsf"),Culonglong,(Ptr{UInt8},Cint),file.rsf,dim)
end

function ucharread(arr::Array{UInt8,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_ucharread,"libdrsf"),Void,(Ptr{UInt8},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function charread(arr::Array{UInt8,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_charread,"libdrsf"),Void,(Ptr{UInt8},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function intread(arr::Array{Int32,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_intread,"libdrsf"),Void,(Ptr{Cint},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function floatread(arr::Array{Float32,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_floatread,"libdrsf"),Void,(Ptr{Cfloat},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function complexread(arr::Array{Complex64,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_complexread,"libdrsf"),Void,(Ptr{Complex64},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function shortread(arr::Array{Int16,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_shortread,"libdrsf"),Void,(Ptr{Cshort},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function ucharwrite(arr::Array{UInt8,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_ucharwrite,"libdrsf"),Void,(Ptr{UInt8},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function charwrite(arr::Array{UInt8,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_charwrite,"libdrsf"),Void,(Ptr{UInt8},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function intwrite(arr::Array{Int32,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_intwrite,"libdrsf"),Void,(Ptr{Cint},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function floatwrite(arr::Array{Float32,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_floatwrite,"libdrsf"),Void,(Ptr{Cfloat},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function complexwrite(arr::Array{Complex64,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_complexwrite,"libdrsf"),Void,(Ptr{Complex64},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function shortwrite(arr::Array{Int16,1},size::Integer,file::File)
    size::Csize_t = size
    ccall((:sf_complexwrite,"libdrsf"),Void,(Ptr{Cshort},Csize_t,Ptr{UInt8}),arr,size,file.rsf)
end

function putint(file::File,name::String,val::Integer)
    val::Cint = val
    ccall((:sf_putint,"libdrsf"),Void,(Ptr{UInt8},Ptr{UInt8},Cint),file.rsf,name,val)
end

function putfloat(file::File,name::String,val::Real)
    val::Cfloat = val
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

The size of `file`, an Int array representing the length of each of its
dimensions.

# Examples

```julia-repl
julia> run(pipeline(`sfspike n1=2 n2=3`, stdout="spike.rsf"))

julia> inp = m8r.input("spike.rsf")

julia> m8r.size(inp)
(2, 3)
```
"""
function size(file::File)
    size = leftsize(file, 0)
    dim = 1
    n = histint(file, string("n", dim))
    s = [n]
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
    rsf_read(file)

Reads RSF `file`, returning its contents and header. `file` may be file handle 
(`m8r.File`), filename (`m8r.File.tag`), or `NTuple{2, Base.PipeEndpoint}`. This
last option is useful for reading pipes from `sf` commands (see last example).
When called with `headers_only` keyword argument set to `true`, does not read
contents, only headers.

# Examples

## Reading file handle
```julia-repl
julia> run(pipeline(`sfspike n1=2 n2=3`, stdout="spike.rsf"))

julia> inp = input("spike.rsf")

julia> dat, n, d, o, l, u = rsf_read(inp)
(Float32[1.0 1.0 1.0; 1.0 1.0 1.0], [2, 3], Float32[0.004, 0.1], Float32[0.0, 0.0], String["Time", "Distance"], String["s", "km"])
```

## Reading file name
```julia-repl
julia> rsf_read("spike.rsf")
(Float32[1.0 1.0 1.0; 1.0 1.0 1.0], [2, 3], Float32[0.004, 0.1], Float32[0.0, 0.0], String["Time", "Distance"], String["s", "km"])
```

## Reading pipe
```julia-repl
julia> sfspike(;n1=1, n2=2) |> rsf_read
(Float32[1.0 1.0 1.0; 1.0 1.0 1.0], [2, 3], Float32[0.004, 0.1], Float32[0.0, 0.0], String["Time", "Distance"], String["s", "km"])
```
"""
function rsf_read(file::File; headers_only::Bool=false)
    types = [
         UInt8, # SF_UCHAR
         UInt8, # SF_CHAR
         Int32, # SF_INT
         Float32, # SF_FLOAT
         Complex64, # SF_COMPLEX
         Int16, # SF_SHORT
         Float64, # SF_DOUBLE
         Clong, # SF_LONG (UNIX: Int, Windows: Int32)
        ]
    n = Int[i for i in size(file)]
    sz = prod(n)

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

    if headers_only
        return n, d, o, l, u
    end

    itypes = gettype(file)
    data = zeros(types[itypes], sz)
    if itypes == 1
        ucharread(data, sz, file)
    elseif itypes == 2
        charread(data, sz, file)
    elseif itypes == 3
        intread(data, sz, file)
    elseif itypes == 4
        floatread(data, sz, file)
    elseif itypes == 5
        complexread(data, sz, file)
    elseif itypes == 6
        shortread(data, sz, file)
    else
        throw("Cannot read long, double")
    end

    data = reshape(data, n...)

    return data, n, d, o, l, u
end

rsf_read(name::String; headers_only::Bool=false) =
    rsf_read(input(name); headers_only=headers_only)

function rsf_read(stdin::NTuple{2, Base.PipeEndpoint}; headers_only::Bool=false)
    rin, win = stdin
    flush(win)
    old_stdin = STDIN
    redirect_stdin(rin)
    data = rsf_read("in"; headers_only=headers_only)
    redirect_stdin(old_stdin)
    return data
end

"""
    rsf_write(file, dat, n, d, o, label, unit)

Write RSF `file`. `file` may be file handle (`m8r.File`), filename
(`m8r.File.tag`) or absent:

    rsf_write(dat, n, d, o, label, unit) -> NTuple{2, Base.PipeEndpoint}

This last option is useful for writing pipes to `sf` commands (see last
example). However, it must be noted that in those cases `rsf_write` can be
omitted.

In all methods, `n`, `d`, `o`, `label`, and `unit` are optional. If given, they
should be of type `AbstractArray`.

Finally, one may write from a pipe to a file with:

    rsf_write(file::String, stdin::NTuple{2, Base.PipeEndpoint})

!!! warning "Writing to file handles"

    Do *not* supply the file as an `m8r.File` type unless you know exactly what
    you are doing. Because of how Madagascar is set up, calling `m8r.output`
    automatically sets the filetype to whatever was read in the previous
    `m8r.input("in")` call. If this function has not yet been called, it
    defaults to `float`. Therefore, it is *impossible* to create a `complex`
    file immediately after reading a `float`-type file. We overcome this
    limitation with other `rsf_write` methods by writing dummy files to dummy
    pipes to "trick" Madagascar into switching file types.

    In addition, do *not* try to write to a file handle which has already been
    written to with `rsf_write`, as `rsf_write` closes the file. Doing so will
    cause a segfault.

# Examples

## Write file by name
```julia-repl
julia> rsf_write("spike.rsf", [1., 2.])

julia> rsf_read("spike.rsf")
(Float32[1.0, 2.0], [2], Float32[1.0], Float32[0.0], String[""], String[""])
```

## Write to pipe
```julia-repl
julia> rsf_write([1. im]) |> sfreal |> rsf_read
(Complex{Float32}[1.0+0.0im 0.0+1.0im], [1, 2], Float32[1.0, 1.0], Float32[0.0, 0.0], String["", ""], String["", ""])
```

## Write from pipe
```julia-repl
julia> sfspike(;n1=1) |> x -> rsf_write("spike.rsf", x)

julia> rsf_read("spike.rsf")
(Float32[1.0, 2.0], [2], Float32[1.0], Float32[0.0], String[""], String[""])
```

## Writing file handle (avoid this!)
```julia-repl
julia> out = output("spike.rsf")

julia> rsf_write(out, [1., 2]) # rsf_write(out, [im, 2]) will not work due to warning

julia> rsf_read(out.tag)
(Float32[1.0, 2.0], [2], Float32[1.0], Float32[0.0], String[""], String[""])
```
"""
function rsf_write(file::File, dat::AbstractArray, n=nothing, d=nothing,
                   o=nothing, l=nothing, u=nothing)
    if n == nothing
        n = Base.size(dat)
    end
    dim = length(n)
    d = d == nothing ? [1 for i in 1:dim] : d
    o = o == nothing ? [0 for i in 1:dim] : o
    l = l == nothing ? ["" for i in 1:dim] : l
    u = u == nothing ? ["" for i in 1:dim] : u
    for i in 1:dim
        typeof(n[i]) <: Integer || throw("All n must be `Integer`")
        typeof(d[i]) <: Real    || throw("All d must be `Real`")
        typeof(o[i]) <: Real    || throw("All o must be `Real`")
        typeof(l[i]) <: String  || throw("All l must be `String`")
        typeof(u[i]) <: String  || throw("All u must be `String`")
        putint(file, "n$i", n[i])
        putfloat(file, "d$i", d[i])
        putfloat(file, "o$i", o[i])
        putstring(file, "label$i", l[i])
        putstring(file, "unit$i", u[i])
    end

    if eltype(dat) <: UInt8
        charwrite(Array{UInt8}(vec(dat)), leftsize(file, 0), file)
    elseif eltype(dat) <: AbstractFloat
        floatwrite(Array{Float32}(vec(dat)), leftsize(file, 0), file)
    elseif eltype(dat) <: Complex
        complexwrite(Array{Complex64}(vec(dat)), leftsize(file, 0), file)
    elseif eltype(dat) <: Int16
        shortwrite(Array{Int16}(vec(dat)), leftsize(file, 0), file)
    elseif eltype(dat) <: Integer
        intwrite(Array{Int32}(vec(dat)), leftsize(file, 0), file)
    else
        throw("Cannot write long, double")
    end
    close(file)
end

function rsf_write(name::String, dat::AbstractArray, n=nothing, d=nothing,
                   o=nothing, l=nothing, u=nothing)
    # Madagascar's output function inherits the type of the previous input.
    # Therefore, in order to have the correct output type, one must create a
    # dummy input of the correct type.
    old_stdin = STDIN
    (rin, win) = redirect_stdin()
    spike = joinpath(RSFROOT, "bin", "sfspike")
    if eltype(dat) <: Int16
        dd = joinpath(RSFROOT, "bin", "sfdd")
        pipe = pipeline(`$spike n1=1`, `$dd type=short`)
    elseif eltype(dat) <: Complex
        rtoc = joinpath(RSFROOT, "bin", "sfrtoc")
        pipe = pipeline(`$spike n1=1`, `$rtoc out=stdout`)
    elseif eltype(dat) <: Integer
        dd = joinpath(RSFROOT, "bin", "sfdd")
        pipe = pipeline(`$spike n1=1`, `$dd type=int`)
    else
        pipe = `$spike n1=1 out=stdout`
    end
    Base.wait(spawn(pipeline(pipe, stdout=win)))
    redirect_stdin(old_stdin)
    rsf_read((rin, win))

    rsf_write(output(name), dat, n, d, o, l, u)
end
rsf_write(file::Union{String, File}, dat::AbstractArray; n=nothing, d=nothing,
          o=nothing, l=nothing, u=nothing) = rsf_write(file, dat, n, d, o, l, u)

function rsf_write(dat::AbstractArray, n=nothing, d=nothing, o=nothing,
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
    rsf_write(name, dat, n, d, o, l, u)

    # 2) Pipe it to dummy sfwindow
    old_stdin = STDIN
    (rin, win) = redirect_stdin()
    progpath = joinpath(RSFROOT, "bin", "sfwindow")
    pipe = `$progpath squeeze=n out=stdout`
    Base.wait(spawn(pipeline(pipe, stdin=name, stdout=win)))
    redirect_stdin(old_stdin)

    # 3) Remove temp
    progpath = joinpath(RSFROOT, "bin", "sfrm")
    run(pipeline(`$progpath $name`))
    dir = dirname(name)
    spawn(pipeline(`rmdir $dir`))
    return rin, win
end
rsf_write(dat::AbstractArray; n=nothing, d=nothing, o=nothing, l=nothing,
          u=nothing) = rsf_write(dat, n, d, o, l, u)

function rsf_write(name::String, stdin::NTuple{2, Base.PipeEndpoint})
    dat = rsf_read(stdin)
    rsf_write(name, dat...)
end

function process_args(;kwargs...)
    args = String[]
    for (key, val) in kwargs
        if typeof(val) <: Tuple
            val = join(["$v" for v in val], ",")
        elseif typeof(val) <: Bool
            val = val ? "y" : "n"
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
            manfile = joinpath(m8r.RSFROOT, "share", "man", "man1",
                               progname*".1")
            if isfile(manfile)
                old_stdout = STDOUT
                (rout, wout) = redirect_stdout()
                run(pipeline(`man $manfile`, stdout=wout, stdin=DevNull,
                             stderr=DevNull))
                Base.close(wout)
                manpage = convert(String, readavailable(rout))
                manpage = replace(manpage, "\n", "\n\t")
                manpage = "\n# RSF Documentation\n"*manpage
                Base.close(rout)
                redirect_stdout(old_stdout)
            else
                manpage = ""
            end

"""
    $progname(input; to_file, kwargs...) -> output

Runs RSF program `$progname` on the data provided by `input`. This may be an
`m8r.File`, a `NTuple{2, Base.PipeEndpoint}` (that is, the output of another
`sf` command) or, an array (and optionally, positional arguments n, d, o, l, u).
If the program requires no input, it may be absent.

If the keyword argument `to_file` is absent, `output` is a 2-tuple of pipes,
that is, `NTuple{2, Base.PipeEndpoint}`. If `to_file` is present, contents are
written to the filename it points to. In this case, `output` is `nothing`.

It is also possible to pass keyword arguments to the `$progname`. See `?m8r`
for examples.

$manpage"""
            function ($F)(stdin::NTuple{2, Base.PipeEndpoint};
                            kwargs...)
                rin, win = stdin
                returnval = (rin, win)
                outstdout = "out=stdout"
                i = find(x->x[1] == :to_file, kwargs)
                if length(i) > 0
                    win = kwargs[i[]][2]
                    deleteat!(kwargs, i[])
                    outstdout = ""
                    returnval = nothing
                end
                args = process_args(;kwargs...)
                progpath = joinpath(RSFROOT, "bin", $S)
                pipe = `$progpath $args $outstdout`
                Base.wait(spawn(pipeline(pipe, stdin=rin, stdout=win)))
                return returnval
            end
        end
        @eval function ($F)(;kwargs...)
                old_stdin = STDIN
                (rin, win) = redirect_stdin()
                returnval = (rin, win)
                outstdout = "out=stdout"
                i = find(x->x[1] == :to_file, kwargs)
                if length(i) > 0
                    win = kwargs[i[]][2]
                    deleteat!(kwargs, i[])
                    outstdout = ""
                    returnval = nothing
                end
                args = process_args(;kwargs...)
                progpath = joinpath(RSFROOT, "bin", $S)
                pipe = `$progpath $args $outstdout`
                Base.wait(spawn(pipeline(pipe, stdout=win)))
                redirect_stdin(old_stdin)
                return returnval
        end
        @eval function ($F)(tag::String; kwargs...)
                old_stdin = STDIN
                (rin, win) = redirect_stdin()
                returnval = (rin, win)
                outstdout = "out=stdout"
                i = find(x->x[1] == :to_file, kwargs)
                if length(i) > 0
                    win = kwargs[i[]][2]
                    deleteat!(kwargs, i[])
                    outstdout = ""
                    returnval = nothing
                end
                args = process_args(;kwargs...)
                progpath = joinpath(RSFROOT, "bin", $S)
                pipe = `$progpath $args $outstdout`
                Base.wait(spawn(pipeline(pipe, stdin=tag, stdout=win)))
                redirect_stdin(old_stdin)
                return returnval
        end
        @eval function ($F)(file::File; kwargs...)
            return $F(file.tag; kwargs...)
        end
        @eval function ($F)(dat::AbstractArray, n=nothing, d=nothing, o=nothing,
            l=nothing, u=nothing; kwargs...)
            return rsf_write(dat, n, d, o, l, u) |> x -> $F(x; kwargs...)
        end
    end
end

end
