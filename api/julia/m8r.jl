module m8r

export sf_init

type File
     tag::ASCIIString
end

function sf_init()
	 ccall((:sf_init,"librsf"),Void,(Int32,Ptr{Ptr{Uint8}}),(length(ARGS),ARGS))
end

end