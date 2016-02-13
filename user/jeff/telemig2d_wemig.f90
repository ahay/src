module telemig2d_wemig
  	use telemig2d_kernel       !! GOGI's 80 degree code
	use fft
	use rsf
	
  	implicit none  
 
  	integer,private :: nx,nz,nw,nh,ih,ix,iz,iw,id,ntaper
  	real,private :: ox,dx,oz,dz,ow,dw,oh,dh,nxtap,pi,eps
  	logical ,private:: source_norm,forward,verbose
  
  	real,allocatable,dimension(:) :: tap,ww,saxmax
  	complex,allocatable,dimension(:,:) :: sax,rax
  	real,allocatable,dimension(:,:,:) :: Timg
  	real,allocatable,dimension(:,:,:,:) :: Tcig

contains
  	!----------------------------------------------------------------
  	!! . . Initialization routine - parameter passing and array allocation
  	subroutine wemig_axes_init(nx_in,ox_in,dx_in,nz_in,oz_in,dz_in,nw_in,ow_in,dw_in,nh_in,oh_in,dh_in,&
  	         source_norm_in,verbose_in,forward_in,ntaper_in,eps_in)
    	integer ::nx_in,nz_in,nw_in,nh_in,ntaper_in
     	real :: ox_in,dx_in,oz_in,dz_in,ow_in,dw_in,oh_in,dh_in,eps_in
        logical :: source_norm_in, verbose_in,forward_in
        
		!! For OMP
    	integer, external :: omp_get_num_threads,omp_get_thread_num
    	integer  :: ix,iz,ih,nth,ith
     	
     	!! . . Distance
     	nx = nx_in; ox = ox_in; dx = dx_in
     
     	!! . . Depth
     	nz = nz_in; oz = oz_in; dz = dz_in
     
       	!! . . Offset
       	nh = nh_in; oh=oh_in; dh=dh_in
     
     	!! . . Frequency
     	pi=acos(-1.)
     	nw=nw_in; dw=2.*pi*dw_in; ow=2.*pi*ow_in
  	 	
  	 	source_norm = source_norm_in
  	 	verbose = verbose_in
  	 	ntaper=ntaper_in
  	 	forward=forward_in
  	 	eps=eps_in

#ifdef _OPENMP
     !$OMP PARALLEL
     nth = omp_get_num_threads()
     !$OMP END PARALLEL
#else
     nth = 1
#endif

    	!! Velocity and fields
    	allocate( tap(nx),sax(nx,nth),rax(nx,nth),ww(nth))
    	allocate( Tcig(nh,nx,nz,nth) )
    	allocate( Timg(nx,nz,nth) )
    	allocate( saxmax(nth) )

    	!$OMP PARALLEL DO PRIVATE(ith,iz,ih,ix) 
   	 	do ith=2,nth
			do iz=1,nz
	        	do ix=1,nx
                	Timg(ix,iz,ith)=0.
              	end do
			end do
		end do
    	!$OMP END PARALLEL DO
 
 		!$OMP PARALLEL DO PRIVATE(ith,iz,ih,ix) 
   	 	do ith=2,nth
			do iz=1,nz
	        	do ix=1,nx
        			do ih=1,nh
                		Tcig(ih,ix,iz,ith)=0.
              		end do
        		end do
			end do
		end do
    	!$OMP END PARALLEL DO
    	
    	!! .. Set up taper
    	tap=1.
    	call rwetaper_init(ntaper) !! Initialize the taper program
  	end subroutine wemig_axes_init
 	
 	!----------------------------------------------------------------
  	!! . . Main wavefield extrapolation and imaging routine
  	subroutine wemig(velS,velR,swf,rwf,img,cig)
    	real    :: rarg,sarg,ss,w,vminS,vminR
   		real    :: velS(:,:), velR(:,:), img(:,:),cig(:,:,:)
    	complex :: swf(:,:),rwf(:,:)

    	!! For OMP
    	integer, external :: omp_get_num_threads,omp_get_thread_num
    	integer  :: nth,ith
    	logical :: dofilt

#ifdef _OPENMP
    	!$OMP PARALLEL
    		nth = omp_get_num_threads()
      !$OMP END PARALLEL
#else
      nth = 1
#endif

		!! . . Source directivity
    	sarg =-1.;    
    	rarg = 1.; if (forward) rarg=-1.

 		!! . . Get the minimum velocity for phase filtering
 		vminS = minval( velS )
 		vminR = minval( velR )
 
 		!! . . Initialize the migration kernel
    	call shot_ker_init(nth, nx, dx, dz, 0., 0., .11)

    	!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iw,id,ix,ih,iz) 
    	Frequency_loop: do iw=1,nw  !! FREQUENCY LOOP       
    
        !! . . Thread ID number
#ifdef _OPENMP
       		id = omp_get_thread_num()+1
#else
         id = 1
#endif
			!! . . Local thread Frequency 
       		ww(id) = real(iw-1.)*dw+ow

       		!! . . Define the Source wavefield 
       		do ix=1,nx
          		sax(ix,id) = swf(ix,iw)
      	 	end do

       		!! . . Define Receiver wavefield 
       		do ix=1,nx
          		rax(ix,id) = rwf(ix,iw)
       		end do
       
      	 	!! . . Taper wavefields
       		call rwetaper(sax,id) 
       		call rwetaper(rax,id) 

       		!! . . Find the max amplitude of the SWF for normalization
       		if (source_norm) saxmax(id) = maxval( real(conjg(sax(:,id))*sax(:,id)) )

       		Extrap_loop: do iz=2,nz  

				!! .. SWF maximum for deconvoluation
				saxmax(id) = maxval (abs ( sax(:,id)))

          		!! . . Step Source and Receiver wavefields
          		call shot_ker_onestep(id, sax, velS(:,iz), ww(id), sarg)
          		call shot_ker_onestep(id, rax, velR(:,iz), ww(id), rarg)

          		!! . . Taper
          		call rwetaper(sax,id) 
          		call rwetaper(rax,id) 

          		!! . . High-angle filter  FFT to kx
          		call fth(.true., .false., sax(:,id) )
          		call fth(.true., .false., rax(:,id) )

          		!! . . Phase correction routine
         		call wemig_phs_correction(id,sax,ww(id),vminS,sarg)
          		call wemig_phs_correction(id,rax,ww(id),vminR,rarg)

          		!! . . FFT to x
          		call fth(.false., .false., sax(:,id) )
          		call fth(.false., .false., rax(:,id) )

          		!! . . Taper
          		call rwetaper(sax,id) 
          		call rwetaper(rax,id) 

          		!! . . Imaging condition
         		if (source_norm) then
            		do ix=1,nx
                		Timg(ix,iz,id) = Timg(ix,iz,id) + real( rax(ix,id)*conjg(sax(ix,id)) ) &
                		/ ( real( sax(ix,id)*conjg(sax(ix,id))) + 0.00000000001 + eps*saxmax(id))
             		end do
          		else 
            		do ix=1,nx
            	   		Timg(ix,iz,id)=Timg(ix,iz,id)+ real( rax(ix,id)*conjg(sax(ix,id)) )
           		  	end do
          		endif
          		
          		!! . . CIGS
          		 if (nh .gt. 0) then
            		do ih=1,nh
            			do ix=1,nx
            	    	   if ( ix-(nh+1)/2+ih .ge. 1 .and. ix-(nh+1)/2+ih .le. nx .and. &
            	            	ix+(nh+1)/2-ih .ge. 1 .and. ix+(nh+1)/2-ih .le. nx  ) then
            		          		Tcig(ih,ix,iz,id) = Tcig(ih,ix,iz,id) +  &
                	      		real( rax(ix-(nh+1)/2+ih,id) * conjg( sax(ix+(nh+1)/2-ih,id)) ) 
                	   		end if
                		end do
             		end do
          		end if

	       	end do Extrap_loop

    	   	if (id .eq. 1 ) write(0,'(a,i5,2f11.6,2f11.6)')' * ST: F%,w,im: ',100*iw/nw,& 
       			ww(id)/2./pi,maxval(Timg(:,:,id))

    	end do Frequency_loop
    	!$OMP END PARALLEL DO

write(0,*) 'IMG OUTPUT',maxval(Timg)

    	!! SUMMATION OVER OUTSIDE IMAGE CONDITION
    	do ith=1,nth
    		!$OMP PARALLEL DO PRIVATE(ix,ih,iz) 
    	   	do iz=1,nz
    	   	      	do ix=1,nx
       		         	img(ix,iz) = img(ix,iz) + Timg(ix,iz,ith)
          		end do
    		end do
       		!$OMP END PARALLEL DO
   		end do

    	!! SUMMATION OVER OUTSIDE CIG 
    	do ith=1,nth
    		!$OMP PARALLEL DO PRIVATE(ix,ih,iz) 
    	   	do iz=1,nz
    	   		do ix=1,nx
    	      		do ih=1,nh
       		         	cig(ih,ix,iz) = cig(ih,ix,iz) + Tcig(ih,ix,iz,ith)
       	     	 	end do
          		end do
    		end do
       		!$OMP END PARALLEL DO
   		end do
   		
   		
   		!$OMP PARALLEL DO PRIVATE(ith,iz,ix) 
    	do ith=1,nth
    		do iz=1,nz
       			do ix=1,nx
       	        	Timg(ix,iz,ith)=0.
       	   		 end do
			end do
		end do
 		!$OMP END PARALLEL DO

   		!$OMP PARALLEL DO PRIVATE(ith,iz,ih,ix) 
    	do ith=1,nth
    		do iz=1,nz
    	   		do ih=1,nh
       		      	do ix=1,nx
       	        	 	Tcig(ih,ix,iz,ith)=0.
       	   		  	end do
    			end do
			end do
		end do
 		!$OMP END PARALLEL DO
 		
    	call shot_ker_release()

	end subroutine wemig

  	!----------------------------------------------------------------
  	!! . . Setup Taper 
  subroutine rwetaper_init(ntap)
    integer :: itap,jtap,ntap
    real :: pi
    pi=acos(-1.)

    if (ntap .gt. 1) then

       do itap=1,ntap
          jtap = abs(ntap-itap)
          tap(   itap  ) = cos(pi/2.* jtap/(ntap-1))
       end do

       do itap=1,ntap
          jtap = abs(ntap-itap)
          tap(nx-itap+1) = cos(pi/2.* jtap/(ntap-1))
       end do

    end if

  end subroutine rwetaper_init
  	!----------------------------------------------------------------
  	!! . . Apply Taper to wavefield
  	subroutine rwetaper(dax,id)
    	complex :: dax(:,:)
    	integer :: ix,id
    	do ix=1,nx
      	 	dax(ix,id) = dax(ix,id) * tap(ix)
    	end do
  	end subroutine rwetaper

end module  telemig2d_wemig
