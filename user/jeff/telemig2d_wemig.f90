! Module for Teleseismic 2D Migration code

!!$  Copyright (C) 2012 The University of Western Australia
!!$  
!!$  This program is free software; you can redistribute it and/or modify
!!$  it under the terms of the GNU General Public License as published by
!!$  the Free Software Foundation; either version 2 of the License, or
!!$  (at your option) any later version.
!!$  
!!$  This program is distributed in the hope that it will be useful,
!!$  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!$  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!$  GNU General Public License for more details.
!!$  
!!$  You should have received a copy of the GNU General Public License
!!$  along with this program; if not, write to the Free Software
!!$  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

module telemig2d_wemig
  	use telemig2d_kernel

  	implicit none  
 
  	integer :: nx,nz,nw,nh,ih,ix,iz,ith,nth,iw
  	real :: ox,dx,oz,dz,ow,dz,oh,dh,nth,nxtap
  	real :: pi,eps
  	logical :: source_norm,forward
  
  	real,allocatable,dimension(:) :: tap,ww,saxmax
  	real,allocatable,dimension(:,:) :: sax,rax
  	real,allocatable,dimension(:,:,:,:) :: Timg

contains
  	!----------------------------------------------------------------
  	!! . . Initialization routine - parameter passing and array allocation
  	subroutine wemig_axes_init(nx_in,ox_in,dz_in,nz_in,oz_in,dz_in,nw_in,ow_in,dw_in,nh_in,oh_in,dh_in)
    	integer ::nx_in,nz_in,nw_in,nh_in
     	real :: ox_in,dx_in,oz_in,dz_in,ow_in,dw_in,oh_in,dh_in

		!! For OMP
    	integer, external :: omp_get_num_threads,omp_get_thread_num
    	integer  :: nth,ith,ix,iz,ih
     	
     	!! . . Distance
     	nx = nx_in;     ox = ox_in;     dx = dx_in
     
     	!! . . Depth
     	nz = nz_in;     oz = oz_in;     dz = dz_in
     
     	!! . . Frequency
     	pi=acos(-1.)
     	nw=nw_in;     dw=2.*pi*dw_in;      ow=2.*pi*ow_in

     	call wemig_param_init()
  	 	
  	 	!$OMP PARALLEL
    		nth = omp_get_num_threads()
    	!$OMP END PARALLEL

    	!! Velocity and fields
    	allocate( tap(nx),sax(nx,nth),rax(nx,nth),ww(nth))
    	allocate( Timg(nx,nh,nz,nth) )
    	allocate( saxmax(nth) )

    	!$OMP PARALLEL DO PRIVATE(ith,iz,ih,ix) 
   	 	do ith=2,nth
			do iz=1,nz
        		do ih=1,nh
	        		do ix=1,nx
                		Timg(ix,ih,iz,ith)=0.
              		end do
        		end do
			end do
		end do
    	!$OMP END PARALLEL DO
    
    	call rwetaper_init() !! Initialize the taper program
  	end subroutine wemig_axes_init
 	
 	!----------------------------------------------------------------
  	!! . . Main wavefield extrapolation and imaging routine
  	subroutine wemig(velS,velR,swf,rwf,img)
    	real    :: rarg,sarg,ss,w,vminS,vminR
   		real    :: velS(:,:), velR(:,:), img(:,:,:)
    	complex :: swf(:,:),rwf(:,:)

    	!! For OMP
    	integer, external :: omp_get_num_threads,omp_get_thread_num
    	integer  :: nth,ith
    	logical :: dofilt

    	!$OMP PARALLEL
    		nth = omp_get_num_threads()
    	!$OMP END PARALLEL

		!! . . Source directivity
    	sarg =-1.;    if (forward) sarg=1.
    	rarg = 1.

 		!! . . Get the minimum velocity for phase filtering
 		vminS = minval( velS )
 		vminR = minval( velR )
 
 		!! . . Initialize the migration kernel
    	call shot_ker_init(nth, nx, dx, dz, 0., 0., .11)

    	!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(iw,id,ix,ih,iz) &
    	SHARED(nw,ow,dw,nz,nh,nx,rarg,sarg,Timg,sax,rax,saxmax,tap)
    	Frequency_loop: do iw=1,nw  !! FREQUENCY LOOP       
    
    		!! . . Thread ID number
       		id = omp_get_thread_num()+1

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
       		if (source_norm) saxmax(id) = maxval( real(conjg(sax(:,id))*sax(:,id)))

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
          		call SAX_Forward(sax,id)
          		call RAX_Forward(rax,id)

          		!! . . Phase correction routine
         		call wemig_phs_correction(id,sax,ww(id),vminS,sarg)
          		call wemig_phs_correction(id,rax,ww(id),vminR,rarg)

          		!! . . FFT to x
          		call SAX_Backward(sax,id)
          		call RAX_Backward(rax,id)

          		!! . . Taper
          		call rwetaper(sax,id) 
          		call rwetaper(rax,id) 

          		!! . . Imaging condition
         		if (nh .eq. 1 .and. source_norm) then
            		do ix=1,nx
                		Timg(ix,1,iz,id) = Timg(ix,1,ig,id) + real( rax(ix,id)*conjg(sax(ix,id)) ) &
                		/ ( real( sax(ix,id)*conjg(sax(ix,id))) + 0.00000000001 + eps*saxmax(id))
             		end do
          		else if (nh1 .eq. 1) then
            		do ix=1,nx
            	   		Timg(ix,1,iz,id)=Timg(ix,1,iz,id)+ real( rax(ix,id)*conjg(sax(ix,id)) )
           		  	end do
          		else
            		do ix=1,nx
            			do ih=1,nh
            	    	   if ( ix-(nh+1)/2+ih .ge. 1 .and. ix-(nh+1)/2+ih .le. nx .and. &
            	            	ix+(nh+1)/2-ih .ge. 1 .and. ix+(nh+1)/2-ih .le. nx  ) then
            		          		Timg(ix,ih,iz,id) = Timg(ix,ih,ix,id) +  &
                	      		real( rax(ix-(nh+1)/2+ih,id) * conjg( sax(ix+(nh+1)/2-ih,id)) ) 
                	   		end if
                		end do
             		end do
          		end if

	       	end do Extrap_loop

    	   	if (verbose .and. id .eq. 1 ) write(0,'(a,i5,i5,i5,2f11.6,2f11.6)')' * ST: F%,is,ns,w,im: ',100*iw/nw,& 
       			ww(id)/2./pi,maxval(Timg(:,:,:,id))

    	end do Frequency_loop
    	!$OMP END PARALLEL DO

    	!! SUMMATION OVER OUTSIDE IMAGE CONDITION
    	do ith=1,nth
    		!$OMP PARALLEL DO PRIVATE(ie3,ih1,ie1) 
    	   	do iz=1,nz
    	      	do ih=1,nh
    	   	      	do ix=1,nx
       		         	img(ix,ih,iz)=img(ix,ih,iz) + Timg(ix,ih,iz,ith)
       	     	 	end do
          		end do
    		end do
       		!$OMP END PARALLEL DO
   		end do

    	!$OMP PARALLEL DO PRIVATE(ith,iz,ih,ix) 
    	do ith=1,nth
    		do i=1,ne3
    	   		do ih=1,nh
       		      	do ix=1,nx
       	        	 	Timg(ix,ih,iz,ith)=0.
       	   		  	end do
    			end do
			end do
		end do
   		!$OMP END PARALLEL DO

    	call shot_ker_release()

	end subroutine wemig

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
