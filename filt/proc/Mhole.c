program Hole
  use sep

  implicit none
  integer                            :: n1, n2, i1, i2
  real, dimension (:,:), allocatable :: pp, maskout
  real                               :: x,y,u,v

  call sep_init ()
  call from_history (n1, n2)
  call sep_close ()

  allocate (pp(n1,n2), maskout(n1,n2) )
  call sep_read (pp)

	maskout = 1.0

  do i1=1,n1 
     do i2=1,n2 
        x = (i1-1.)/n1 - .5
        y = (i2-1.)/n2 - .3
        u =  x+y
        v = (x-y) / 2.
        if( u**2 + v**2 < .15 ) then
					pp(i1,i2) = 0.
					maskout(i1,i2) = 0.
				end if
      end do
  end do

	call to_history( "n1", n1, "maskout")
	call to_history( "n2", n2, "maskout")
	call sep_write( maskout, "maskout")

  call sep_write (pp)
  deallocate (pp)
  call exit (0)
end program Hole

