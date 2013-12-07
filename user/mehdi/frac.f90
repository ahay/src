! Module to implement backus averaging from sfbackus
module frac
  use rsf

  implicit none
  !integer, private :: navg,nz

  contains

!############################################################################
!####################      Compute C for VTI media  #########################

    subroutine VTI_C(vp_up,vs_up,rho_up,epsilon2,gamma2,delta,cvti)
    real, intent(in):: vp_up,vs_up,rho_up,epsilon2,gamma2,delta
    real, allocatable, intent(out):: cvti(:,:)
    
    allocate (cvti(6,6))  
    cvti(1,1)=rho_up*vp_up**2.;
    cvti(2,2)=rho_up*vp_up**2.;
    cvti(3,3)=cvti(1,1)/(2.*epsilon2+1.);
    cvti(4,4)=rho_up*vs_up**2.;
    cvti(5,5)=rho_up*vs_up**2.;
    cvti(6,6)=gamma2*2.*cvti(4,4)+cvti(4,4);
    cvti(1,2)=cvti(1,1)-2.*cvti(6,6);
    cvti(2,1)=cvti(1,2);
    cvti(1,3)=sqrt(delta*2.*cvti(3,3)*(cvti(3,3)-cvti(4,4))+(cvti(3,3)-&
	cvti(4,4))**2.)-cvti(4,4);
    cvti(2,3)=cvti(1,3);
    cvti(3,2)=cvti(1,3);
    cvti(3,1)=cvti(1,3);
    end subroutine VTI_C
!############################################################################
!###########     Compute generalized anisotropy parameters  #################


    subroutine Aniso_parameters(Dex,Dey,Dez,De15,De16,De24,De26,&
	De34,De35,De45,Ddx,&
    Ddy,Ddz,Dgx,Dgy,Dgz,Dchix,Dchiy,Dchiz,cu,cl,Vp1f,Vs1f,Vp2f,Vs2f,rho1f,rho2f)
    
    
    real  :: exd,exu, eyd,eyu, ezd,ezu,e15d,e15u,e16d,e16u 
    real  ::   e24d,e24u,e26d,e26u,e34d,e34u,e35d,e35u 
    real  ::   e45d,e45u,dxd,dxu,dyd,dyu,dzd,dzu,gxd,gxu,gyd,gyu 
    real  :: 	  gzd,gzu,chixd,chixu,chiyd,chiyu,chizd,chizu
    !real  :: 	  rho1ff, rho2ff, Vp1ff, Vp2ff,Vs1ff, Vs2ff
    real, intent(in)::Vp1f,Vs1f,Vp2f,Vs2f,rho1f,rho2f
    real, allocatable, intent(in):: cu(:,:),cl(:,:)
    real, allocatable :: cn(:,:),cn_lower(:,:)
    real, intent(out):: Dex, Dey, Dez,De15,De16
    real, intent(out)::   De24,De26,De34,De35,De45
    real, intent(out)::   Ddx,Ddy,Ddz,Dgx,Dgy
    real, intent(out):: 	  Dgz,Dchix,Dchiy,Dchiz



!    allocate (cu(6,6))
!    allocate (cl(6,6))
    allocate (cn(6,6))
    allocate (cn_lower(6,6))

    
    
    cn=cu/rho1f;
    exu=(cn(1,1)-Vp1f**2.)/(2*Vp1f**2.);
    eyu=(cn(2,2)-Vp1f**2.)/(2*Vp1f**2.);
    ezu=(cn(3,3)-Vp1f**2.)/(2*Vp1f**2.);
    e15u=(cn(1,5))/(Vp1f**2.);
    e16u=(cn(1,6))/(Vp1f**2.);
    e24u=(cn(2,4))/(Vp1f**2.);
    e26u=(cn(2,6))/(Vp1f**2.);
    e34u=(cn(3,4))/(Vp1f**2.);
    e35u=(cn(3,5))/(Vp1f**2.);
    e45u=(cn(4,5))/(Vs1f**2.);
    
    
    dxu=(cn(1,3)+2.*cn(5,5)-Vp1f**2.)/(Vp1f**2.);
    dyu=(cn(2,3)+2.*cn(4,4)-Vp1f**2.)/(Vp1f**2.);
    dzu=(cn(1,2)+2.*cn(6,6)-Vp1f**2.)/(Vp1f**2.);
    
    gxu=(cn(5,5)-Vs1f**2.)/(2*Vs1f**2.);
    gyu=(cn(4,4)-Vs1f**2.)/(2*Vs1f**2.);
    gzu=(cn(4,4)-Vs1f**2.)/(2*Vs1f**2.);
    
    
    chixu=(cn(1,4)+2.*cn(5,6))/(Vp1f**2.);
    chiyu=(cn(2,5)+2.*cn(4,6))/(Vp1f**2.);
    chizu=(cn(3,6)+2.*cn(4,5))/(Vp1f**2.);
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cn_lower=cl/rho2f
    
    exd=(cn_lower(1,1)-Vp2f**2.)/(2*Vp2f**2.)
    eyd=(cn_lower(2,2)-Vp2f**2.)/(2*Vp2f**2.)
    ezd=(cn_lower(3,3)-Vp2f**2.)/(2*Vp2f**2.)
    e15d=(cn_lower(1,5))/(Vp2f**2.)
    e16d=(cn_lower(1,6))/(Vp2f**2.)
    e24d=(cn_lower(2,4))/(Vp2f**2.)
    e26d=(cn_lower(2,6))/(Vp2f**2.)
    e34d=(cn_lower(3,4))/(Vp2f**2.)
    e35d=(cn_lower(3,5))/(Vp2f**2.)
    e45d=(cn_lower(4,5))/(Vs2f**2.)
    
    
    dxd=(cn_lower(1,3)+2.*cn_lower(5,5)-Vp2f**2.)/(Vp2f**2.);
    dyd=(cn_lower(2,3)+2.*cn_lower(4,4)-Vp2f**2.)/(Vp2f**2.);
    dzd=(cn_lower(1,2)+2.*cn_lower(6,6)-Vp2f**2.)/(Vp2f**2.);
    
    gxd=(cn_lower(5,5)-Vs2f**2.)/(2*Vs2f**2.);
    gyd=(cn_lower(4,4)-Vs2f**2.)/(2*Vs2f**2.);
    gzd=(cn_lower(4,4)-Vs2f**2.)/(2*Vs2f**2.);
    
    
    chixd=(cn_lower(1,4)+2.*cn_lower(5,6))/(Vp2f**2.);
    chiyd=(cn_lower(2,5)+2.*cn_lower(4,6))/(Vp2f**2.);
    chizd=(cn_lower(3,6)+2.*cn_lower(4,5))/(Vp2f**2.);
    
    
    
    Dex=exd-exu;
    Dey=eyd-eyu;
    Dez=ezd-ezu;
    De15=e15d-e15u;
    De16=e16d-e16u;
    De24=e24d-e24u;
    De26=e26d-e26u;
    De34=e34d-e34u;
    De35=e35d-e35u;
    De45=e45d-e45u;
    
    Ddx=dxd-dxu;
    Ddy=dyd-dyu;
    Ddz=dzd-dzu;
    
    Dgx=gxd-gxu;
    Dgy=gyd-gyu;
    Dgz=gzd-gzu;
    
    Dchix=chixd-chixu;
    Dchiy=chiyd-chiyu;
    Dchiz=chizd-chizu;

    end subroutine Aniso_parameters

!############################################################################
!###############################     CG  ####################################
  !! Initialization subroutine
    subroutine cg(nun,const,FFmat,tolerance,itermax,It,w)
      integer     :: It, itermax
      integer                  ::  i11, j11
      real, intent(in)         :: tolerance
      integer, intent(in)         :: nun
      real                         :: alpha, r, beta
      real, allocatable, intent(in) ::   const(:), FFmat(:,:)
      real, allocatable, intent(out) ::  w(:)
      real, allocatable     ::   pp(:), g(:), AA(:,:), B(:), z(:), d(:)

  
 	 allocate (g (nun))
  	 allocate (pp (nun))
 	 allocate (B (nun))
 	 allocate (AA (nun,nun))
 	 allocate (z (nun))
 	 allocate (d (nun))
         allocate (w (nun)); w=0.;
	B=const
	g=-B
	pp=g
	AA=FFmat

	do i11=1,It
    	alpha=dot_product(g,g)/dot_product(pp,matmul(AA,pp))
    	z=dot_product(g,g)
    	w=w-alpha*pp
    	g=g-alpha*matmul(AA,pp)
    	d=dot_product(g,g)
    	beta=sum(d)/sum(z)

   	pp=g+beta*pp
	end do
	
    
        g=matmul(AA,w)-const
        pp=g
        
        do j11=1,itermax
        alpha=dot_product(g,g)/dot_product(pp,matmul(AA,pp))
        z=dot_product(g,g)
        w=w-alpha*pp
        g=g-alpha*matmul(AA,pp)
        d=dot_product(g,g)
        beta=sum(d)/sum(z)
        pp=g+beta*pp
        end do
    
    end subroutine cg
!###############################     CG  ####################################
!############################################################################




!############################################################################
!###############################     Random Vector  #########################

!     Return a random vector of size cc1 with mean=0, max=1 and min=-1
        subroutine RandomV(cc1,randomnoise)

	INTEGER ::  Kseed
        INTEGER, intent(in)  :: cc1
        integer, dimension(8) ::   values
        integer, dimension(1) ::   seed
        real, allocatable ::   harvest(:)
        real, allocatable, intent(out) ::   randomnoise(:)
	allocate (harvest(cc1))
	allocate (randomnoise(cc1))
	call date_and_time(values=values)
	kseed=1
	call random_seed(size=kseed)

	seed= values(8)
	call random_seed(put=seed(1:kseed))
	call random_number(harvest)
        randomnoise=(harvest-.5)*2.

        end subroutine RandomV
!###############################     Random Vector  #########################
!############################################################################




end module frac
