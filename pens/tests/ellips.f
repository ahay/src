c     draw some ellipses  
c
	parameter ( ne = 6 )
	parameter ( nx = 100 , ny = 100 )
c
	real a(ne) , b(ne)
	real xco(2*nx) , yco(2*ny)
	integer col(ne)

	data a/1.,2.,3.,1.,1.,1./
	data b/1.,1.,1.,2.,3.,4./
	data col/6,5,4,3,2,1/
	
	call vp_init()
	call vp_erase()
c      .75 = SCREEN_RATIO from vplot.h
	call vp_stretch(-10.,-10.*.75,10.,10.*.75)
	call vp_umove(-10.,-10.*.75)
	call vp_udraw(10.,-10.*.75)
	call vp_udraw(10.,10.*.75)
	call vp_udraw(-10.,10.*.75)
	call vp_udraw(-10.,-10.*.75)

	do 1 ie =1,ne
	call vp_color( col(ie) )

	xmin = -a(ie)
	dx = 2*abs(xmin)/(nx-1)
c
	call vp_umove(xmin,0.)
	do 10 i = 1,nx
	x = ( xmin + (i-1)*dx )
	y = b(ie)*sqrt ( 1. - x*x/( a(ie)**2 ) )
	xco(i)=x
	yco(i)=y
	call vp_udraw(x,y)
10     continue

	call vp_umove(xmin,0.)
	do 20 i = 1,nx
	x = ( xmin + (i-1)*dx )
	y = - b(ie)*sqrt ( 1. - x*x/( a(ie)**2 ) )
	xco(nx+i)=x
	yco(nx+i)=y
	call vp_udraw(x,y)
20     continue

	 if ( ie .eq. 1 ) call vp_uarea(xco,yco, 2*nx, 1, 4, 4)

1	continue

	stop
	end
