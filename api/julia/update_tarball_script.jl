# binaries = read(`ls $(ENV["RSFROOT"])/bin`, String)
# shared_objects = read(`ls $(ENV["RSFROOT"])/lib`, String)
# u = [
#     replace(e, ".so" => s -> "")
#     for e in split(shared_objects, "\n")
#     if occursin(".so", e)
# ]

# v = ["LibraryProduct(\"$e\", :$e)" for e in u]
# w = join(v, "\n")
# final = replace(w, "++)" => s -> "pp)")
# print(final)

function get_target(; dir::String, root::String=ENV["RSFROOT"], ignore::Array{String}=[])
    branches = Dict("bin" => ("ExecutableProduct", ""), "lib" => ("LibraryProduct", ".so"))
    type, ext = branches[dir]
    files = read(`ls $root/$dir`, String)
    u = [
        replace(e, ext => s -> "")
        for e in split(files, "\n")
        if (e != "") && (ext == "" || occursin(ext, e))
    ]
    u = [e for e in u if !(e in ignore)]

    v = ["$type(\"$e\", :$(replace(e, "-" => "")))" for e in u]
    w = join(v, ",\n    ")
    final = replace(w, "++)" => s -> "pp)")
    return final
end

function main(output_file::String; lib_unsupported::AbstractArray{String}=[], bin_unsupported::AbstractArray{String}=[])
    lib_code = get_target(dir="lib", ignore=lib_unsupported)
    bin_code = get_target(dir="bin", ignore=bin_unsupported)

    final_build_string = "[\n    $lib_code,\n    $bin_code\n]"

    url = "https://github.com/ahay/src.git"
    latest_remote_hash = read(`git ls-remote $url HEAD`, String)
    latest_remote_hash = split(latest_remote_hash)[1]

    curr_dir = dirname(@__FILE__)
    base_file_contents = read("$curr_dir/tarball_base.jl", String)


    final_file_contents = replace(
        base_file_contents,
        "\"!?MADAGASCAR_HASH!?\"" => "\"$latest_remote_hash\"",
        "\"!?DEPENDENT_EXECUTABLES!?\"" => final_build_string)

    open(output_file, "w") do f
        write(f, final_file_contents)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    lib_unsupported::Array{String} = []
    bin_unsupported::Array{String} = ["sfpen", "sfvr3d"]
    main("build_tarballs.jl", lib_unsupported=lib_unsupported, bin_unsupported=bin_unsupported)
end

