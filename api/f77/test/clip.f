	program Clipit
	implicit none
	integer n1, n2, i1, i2
	integer*8 in, out
	integer*8 sf_input, sf_output, sf_leftsize, sf_gettype
	logical sf_getfloat, sf_histint
	real clip, trace(1000)

	call sf_init()
	in = sf_input("in")
	out = sf_output("out")

	if (3 .ne. sf_gettype(in)) 
     &  call sf_error("Need float input")

	if (.not. sf_histint(in,"n1",n1)) then
	   call sf_error("No n1= in input")
	else if (n1 > 1000) then
	   call sf_error("n1 is too long")
	end if
	n2 = sf_leftsize(in,1)

	if (.not. sf_getfloat("clip",clip)) 
     &  call sf_error("Need clip=")

	do 10 i2=1, n2
	   call sf_floatread(trace,n1,in)

	   do 20 i1=1, n1
	      if (trace(i1) >  clip) then
		 trace(i1)=clip
	      else if (trace(i1) < -clip) then
		 trace(i1)=-clip
	      end if
 20	   continue

	   call sf_floatwrite(trace,n1,out)
 10	continue

	stop
	end
