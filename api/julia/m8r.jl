module m8r

export sf_init

immutable File
     tag::ASCIIString
     rsf::Ptr{Uint8}
end

function init()
	 argv = [basename(Base.source_path())]
	 append!(argv,ARGS)
	 ccall((:sf_init,"libdrsf"),Void,(Int32,Ptr{Ptr{Uint8}}),length(argv),argv)
end

function input(tag::ASCIIString)
	 rsf = ccall((:sf_input,"libdrsf"),Ptr{Uint8},(Ptr{Uint8},),tag)
	 File(tag,rsf)
end

function histint(file::File,name::ASCIIString)
	 val = Cint[0]
	 ccall((:sf_histint,"libdrsf"),Bool,(Ptr{Uint8},Ptr{Uint8},Ptr{Cint}),file.rsf,name,val)
	 return val[]
end

function leftsize(file::File,dim::Int)
	 ccall((:sf_leftsize,"libdrsf"),Culonglong,(Ptr{Uint8},Cint),file.rsf,dim)
end

function getfloat(name::ASCIIString)
	 val = Cfloat[0]
	 ccall((:sf_getfloat,"libdrsf"),Bool,(Ptr{Uint8},Ptr{Cfloat}),name,val)
	 return val[]
end

function floatread(arr::Array,size::Int,file::File)
	 ccall((:sf_floatread,"libdrsf"),Void,(Ptr{CFloat},Csize_t,Ptr{Uint8}),arr,size,file.rsf)
end

end