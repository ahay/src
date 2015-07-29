c==============================================================
c	FORTRAN routine displaying polymarkers
c==============================================================
c
	dimension xpoly(20),ypoly(20)
	integer marktype(8),color(8)
	integer   out
	data marktype /2, 3, 4, 5, 20, 21, 22, 23/
	data color /0, 1, 2, 3, 4, 5, 6, 7/

	call vp_init()
c
c  choose standard style for vplot
c
	call vp_style(STANDARD)
c
c  set user window
c
	xmin = 0.
	ymin = 0.
	xmax = 100.
	ymax = 100.
	call vp_stretch(xmin,ymin,xmax,ymax)


	ncol   	= 4
	nrow	= 8
	dcol	= (xmax-xmin)/(ncol+1)
	drow	= (ymax-ymin)/(nrow+1)
	isiz	= 2
	icolor  = 0
c
c  compute positions of markers
c
	do irow = 1,nrow
	   ix = 0
	   iy = 0
	   ypos = ymin + drow*irow
c
c  set marker type,color and size
c
	   ityp = ityp+1
	   mtype= marktype(ityp)
	   icolor = icolor+1
	   isiz = isiz+2
	   do icol =1,ncol
		ix = ix+1
		iy = iy+1
		xpoly(ix) = xmin + dcol*icol
		ypoly(iy) = ypos
	   end do
c
c  display markers
c
	   call vp_color(color(icolor))
	   call vp_pmark(ncol,mtype,isiz,xpoly,ypoly)
	end do
c
c  end of session
c
	stop
	end


