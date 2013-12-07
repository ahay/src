!  Forward modeling and inversion of converted wave amplitudes for fracture parameters
!
!
!
!
!!$  Copyright (C) 2013 UT Austin
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


program converted
  use rsf
  use frac 
  implicit none
  

!#################################################################################################################################
!#################################################################################################################################
!#################################################################################################################################


  
  type (file)   ::    aw, w1, refl,NoisyRps, Winvpred
  type(axa)  :: at,az   ! cube axes

  integer ::  n1, n2, itmax, of, azi, jjj, hhh, cc, p, t, azimuth, offset, I, nu, nr, theta,phi, kseed
  real, parameter :: pi = 3.1415927
  integer, allocatable :: values(:)
  real, allocatable :: F11(:,:), F12(:,:) , F22(:,:), F1111(:,:), F1112(:,:), F1122(:,:), F1222(:,:), F2222(:,:)
  real, allocatable :: rV2_cor(:,:), rV20(:,:), a1(:), b1(:), c1(:), d1(:), e1(:), Rpp0(:), Rpp(:,:), Rpp_FRAC(:,:)
  real, allocatable :: Rpp_aniso(:,:), albet(:), cu(:,:), c(:,:)
  real, allocatable :: f1(:), g1(:), h1(:), q1(:), r1(:), s1(:), t1(:), v1(:), z1(:), w1d(:), w2d(:,:)
  real, allocatable :: F(:,:),  FF(:,:),F1D(:), cons(:), weight1(:), FFT(:,:), noisetest(:)
  real ::  BT, BN, Aa1, AA2, Bb1, Bb2,SWsplitting
  real ::  rho1, vp1, vp2, vs1, ratio, rho2, vs2, epsilon1, delta1, gamma1, eps, gam, divide, deltafix, del, unce
  real ::  devVp2, devVs2, devrho2,dev_eps, dev_del, dev_gam, vvvp, vvvs, rrrho2, snr, eta, frac_area, coef,C11
  real :: C12, C13, C22, C23, C33, C44, C45, C55, C66, C16, C26, C36
  real :: phi1,phi2, Dphi, ct, tet,  x1, x2, x3, tolerancer, m, L
  real :: alpha11, alpha12, alpha22, beta1111, beta1112, beta1122,beta1222,beta2222 
  integer, allocatable ::  seed(:)
  real :: numb
  real, allocatable::noise(:),noise1(:,:), noise2(:,:),Rpp0_background(:,:),Rpp_aniso1(:),Rpp_FRAC1(:),Rpp_Minus_iso_background(:,:)
  real, allocatable :: Rpp_noisy(:,:), Rpp_Noisy_Minus_background(:,:), Rpp_remainder(:,:),Rpp_remainder1(:)
  real :: Dexm, Deym, Dezm,De15m,De16m
  real ::   De24m,De26m,De34m,De35m,De45m
  real ::   Ddxm,Ddym,Ddzm,Dgxm,Dgym
  real :: 	  Dgzm,Dchixm,Dchiym,Dchizm, L1, m1,phi_s1_actual,phi_s1_predicted222,ratio1, vpave,vsave


   call sf_init()            ! initialize RSF
   aw = rsf_output("actual_weights") ! common_offset is the output file name for common offset section
   w1 = rsf_output("weights") ! common_offset is the output file name for common offset section
   Refl = rsf_output("Rps") ! common_offset is the output file name for common offset section
   NoisyRps = rsf_output("NoisyRps") ! common_offset is the output file name for common offset section
   Winvpred = rsf_output("Wip") ! common_offset is the output file name for common offset section


   call from_par("phi",phi)   ! command-line parameter
   call from_par("theta",theta)   ! command-line parameter

  call to_par(aw,"n1",8)     ! write the sample number for the 2nd axis for the output file
  call to_par(aw,"o1",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(aw,"d1",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file

  call to_par(w1,"n1",8)     ! write the sample number for the 2nd axis for the output file
  call to_par(w1,"o1",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(w1,"d1",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file


  call to_par(refl,"n1",theta)     ! write the sample number for the 2nd axis for the output file
  call to_par(Refl,"o1",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(Refl,"d1",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file
  call to_par(refl,"n2",phi)     ! write the sample number for the 2nd axis for the output file
  call to_par(refl,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(refl,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file

  call to_par(NoisyRps,"n1",theta)     ! write the sample number for the 2nd axis for the output file
  call to_par(NoisyRps,"o1",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(NoisyRps,"d1",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file
  call to_par(NoisyRps,"n2",phi)     ! write the sample number for the 2nd axis for the output file
  call to_par(NoisyRps,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(NoisyRps,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file
  

  call to_par(Winvpred,"n1",8)     ! write the sample number for the 2nd axis for the output file
  call to_par(Winvpred,"o1",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(Winvpred,"d1",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file
  call to_par(Winvpred,"n2",2)     ! write the sample number for the 2nd axis for the output file
  call to_par(Winvpred,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(Winvpred,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file

  rho1=2530.;      !%density of upper medium
  Vp1=4509.;       !%Vp of upper medium
  Vs1=2855.;       !%Vs of upper medium
  
  
  ratio1=1;
  rho2=2460.*ratio1;      !%density of lower medium
  Vp2=4161.*ratio1;       !%Vp of upper medium
  Vs2=2687.*ratio1;      !%Vs of upper medium
  
  epsilon1=.1;
  delta1=0.1;
  gamma1=.10;
  !%%%%%%% anisotropic parameters for lower  %%%%%%%%%%%%
  eps=0.290;
  del=0.17;
  gam=.1;
  divide=100;
   deltafix=del
  unce=0.15;    !% uncertainty in Monte Carlo simulation default: 15%
  
  !%%%%%%%%%%%%%%%%%%%%%%%%   Generating random background properties %%%%%%%%%%%%%%%%%%%%
  
!  devVp2=ceiling(Vp2*unce);
 ! devVs2=ceil(Vs2*unce);
 ! devrho2=ceil(rho2*unce);
  
!  dev_eps=ceil(100*eps*unce);
!  dev_del=ceil(100*del*unce);
!  dev_gam=ceil(100*gam*unce);
  
!  vvvp=Vp2;
!  vvvs=Vs2;
!  rrrho2=rho2;
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  SNR=2.;             !%Signal to Noise ratio of reflectivity data
  eta=.3;             !%fracture density
  ratio=.75;           !% Bn/Bt ratio
  BT=9.*10.**(-8.)         !%tangential compliance of fracctures
  BN=ratio*BT;        !%normal comliance of fracctures
  frac_area=.00016;       !%10cm x 1 mm   set to give about 10% shear wave splitting
  coef=2.;             !%a constant to adjust values of alpha and betas in forward modeling
  
  
  phi1=50.;     ! %azimuth of the first fracture set
  phi2=-30.;      !%azimuth of the second fracture set
  Dphi=phi2-phi1;
  
  I=20;               !%number of iterations to get good initial estimates for Conjugate Gradient
  itmax=10;        ! %number of iterations for Conjugate Gradient
  !% del=3;
  theta=35;           !%range of offset in data 0-theta
  phi=360;             !%range of azimuth in data 0-phi
  
  

  
  VpAve=(Vp1+Vp2)/2.;
  VsAve=(Vs1+Vs2)/2.;
  
  
  m1=Vs1**2.*rho1;
  L1=rho1*Vp1**2.-2.*m1;
  
  allocate(c(6,6))
  allocate(cu(6,6))
    call  VTI_C(vp1,vs1,rho1,epsilon1,gamma1,delta1,cu)
    call  VTI_C(vp2,vs2,rho2,eps,gam,del,c)
  
  
  
   C11=c(1,1);
   C22=c(2,2);
   C33=c(3,3);
   C12=c(1,2);
   C13=c(1,3);
   C23=c(2,3);
   C44=c(4,4);
   C55=c(5,5);
   C66=c(6,6);
  
   m=Vs2**2.*rho2;
   L=rho2*Vp2**2.-2.*m;
  
  
  
  !#######################################################################################
  !###############################   alpha-beta      #####################################
  

  Aa1=eta*BT*frac_area;
  Aa2=(1.-eta)*BT*frac_area;
  
  
  Bb1=eta*(BN-BT)*frac_area;
  Bb2=(1.-eta)*(BN-BT)*frac_area;
  
  
  phi1=phi1*pi/180.;
  phi2=phi2*pi/180.;
  
  
  
  alpha11=coef*(Aa1*sin(phi1)*sin(phi1)+Aa2*sin(phi2)*sin(phi2));
  alpha12=coef*(Aa1*cos(phi1)*sin(phi1)+Aa2*cos(phi2)*sin(phi2));
  alpha22=coef*(Aa1*cos(phi1)*cos(phi1)+Aa2*cos(phi2)*cos(phi2));
  
  
  beta1111=coef*(Bb1*sin(phi1)**4+Bb2*sin(phi2)**4);
  beta1112=coef*(Bb1*cos(phi1)*sin(phi1)**3+Bb2*cos(phi2)*sin(phi2)**3);
  beta1122=coef*(Bb1*cos(phi1)**2*sin(phi1)**2+Bb2*cos(phi2)**2*sin(phi2)**2);
  beta1222=coef*(Bb1*cos(phi1)**3*sin(phi1)+Bb2*cos(phi2)**3*sin(phi2));
  beta2222=coef*(Bb1*cos(phi1)**4+Bb2*cos(phi2)**4);
  
  
  SWsplitting=((alpha22-alpha11)*C55**2)/(2*C55-(alpha11+alpha22)*C55**2)
  
  
  !###############################   alpha-beta      #####################################
  !#######################################################################################
  

  allocate (F11(theta,phi))
  allocate (F12(theta,phi))
  allocate (F22(theta,phi))
  allocate (F1111(theta,phi))
  allocate (F1112(theta,phi))
  allocate (F1122(theta,phi))
  allocate (F1222(theta,phi))
  allocate (F2222(theta,phi))


  allocate (Rpp0_background(theta,phi))
  allocate (Rpp_Minus_iso_background(theta,phi))
  allocate (Rpp_remainder(theta,phi))
  allocate (Rpp_Noisy_Minus_background(theta,phi))
  allocate (Rpp_noisy(theta,phi))
  allocate (Rpp0(theta))
  allocate (Rpp(theta,phi))
  allocate (Rpp_FRAC(theta,phi))
  allocate (Rpp_aniso(theta,phi))


    call Aniso_parameters(Dexm,Deym,Dezm,De15m,De16m,De24m,De26m,&
	De34m,De35m,De45m,Ddxm,&
    Ddym,Ddzm,Dgxm,Dgym,Dgzm,Dchixm,Dchiym,Dchizm,cu,c,Vp1,Vs1,Vp2,Vs2,rho1,rho2)




  Do t=1,theta
      
      Rpp0(t) = ((rho2*Vp2 - rho1*Vp1)/(rho2*Vp2 +rho1*Vp1)) + (((Vp2 - Vp1)/(Vp2 + Vp1))-&
          4.*((Vs2 + Vs1)**2)*(m - m1)/(((Vp2 + Vp1)**2)*(m + m1)))*(sin(t*pi/180)**2.) + ((Vp2 - &
          Vp1)/(Vp2 + Vp1))*(sin (t*pi/180.)**2.)*(tan(t*pi/180)**2.);
      Do p=1,phi
          
          
          
          
          F11(t,p)=-(C13**2./(4.*(L+2.*m)))+&
              (1./2.)*(cos((p*pi)/180.)**2.*(-((C11*C13+2.*C55**2.)/(L+2.*m))+&
              (4.*C55**2.*VsAve**2.)/(m*VpAve**2.))-&
              (C12*C13*sin((p*pi)/180.)**2.)/(L+2.*m)+C13**2./(2.*(L+2.*m)))*sin((pi*t)/180.)**2.+&
              (1./2.)*(-((C11**2.*cos((p*pi)/180.)**4.)/(2.*(L+2.*m)))-&
              ((C11*C12+2.*C66**2.)*cos((p*pi)/180.)**2.*sin((p*pi)/180.)**2.)/&
              (L+2.*m)-(C12**2.*sin((p*pi)/180.)**4.)/(2.*(L+2.*m)))*&
              sin((pi*t)/180.)**2.*tan((pi*t)/180.)**2.;
          
          F12(t,p)=cos((p*pi)/180.)*((-2.*C55**2.-2.*C13*C66)/(L+2.*m)+&
              (4.*C55**2.*VsAve**2.)/(m*VpAve**2.))*sin((p*pi)/180.)*&
              sin((pi*t)/180.)**2.+cos((p*pi)/180.)*sin((p*pi)/180.)*&
              (-(((C11+C12)*C66*cos((p*pi)/180.)**2.)/(L+2.*m))-&
              ((C11+C12)*C66*sin((p*pi)/180.)**2.)/(L+2.*m))*&
              sin((pi*t)/180.)**2.*&
              tan((pi*t)/180.)**2.;
          
          F22(t,p)=-(C13**2./(4.*(L+2.*m)))+&
              (1./2.)*(-((C12*C13*cos((p*pi)/180.)**2.)/(L+2.*m))+&
              (((-C11)*C13-2.*C55+2.*(C55-C55**2.))/(L+2.*m)+&
              (4.*C55**2.*VsAve**2.)/(m*VpAve**2.))*&
              sin((p*pi)/180.)**2.+C13**2./(2.*(L+2.*m)))*&
              sin((pi*t)/180.)**2.+&
              (1./2.)*(-((C12**2.*cos((p*pi)/180.)**4.)/(2.*(L+2.*m)))-&
              ((C11*C12+2.*C66**2.)*cos((p*pi)/180.)**2.*sin((p*pi)/180.)**2.)/&
              (L+2.*m)-(C11**2.*sin((p*pi)/180.)**4.)/(2.*(L+2.*m)))*&
              sin((pi*t)/180.)**2.*tan((pi*t)/180.)**2.;
          
          F1111(t,p)=-(C13**2./(4.*(L+2.*m)))+&
              (1./2.)*(-((C11*C13*cos((p*pi)/180.)**2.)/(L+2.*m))-&
              (C12*C13*sin((p*pi)/180.)**2.)/(L+2.*m)+C13**2./(2.*(L+2.*m)))*sin((pi*t)/180.)**2.+&
              (1./2.)*(-((C11**2.*cos((p*pi)/180.)**4.)/(2.*(L+2.*m)))-&
              (C11*C12*cos((p*pi)/180.)**2.*sin((p*pi)/180.)**2.)/(L+2.*m)-&
              (C12**2.*sin((p*pi)/180.)**4.)/(2.*(L+2.*m)))*sin((pi*t)/180.)**2.*&
              tan((pi*t)/180.)**2.;
          
          F1112(t,p)=-((2.*C13*C66*cos((p*pi)/180.)*sin((p*pi)/180.)*&
              sin((pi*t)/180.)**2.)/&
              (L+2.*m))+cos((p*pi)/180.)*sin((p*pi)/180.)*&
              (-((2.*C11*C66*cos((p*pi)/180.)**2.)/(L+2.*m))-&
              (2.*C12*C66*sin((p*pi)/180.)**2.)/(L+2.*m))*sin((pi*t)/180.)**2.*&
              tan((pi*t)/180.)**2.;
          
          
          
          F1122(t,p)=-(C13**2./(2.*(L+2.*m)))+&
              (1./2.)*(-(((C11+C12)*C13*cos((p*pi)/180.)**2.)/(L+2.*m))+&
              (((-C11)*C13-C12*C13)*sin((p*pi)/180.)**2.)/(L+2.*m)+C13**2./(L+2.*m))*&
              sin((pi*t)/180.)**2.+&
              (1./2.)*(-((C11*C12*cos((p*pi)/180.)**4.)/(L+2.*m))-&
              ((C11**2.+C12**2.+8.*C66**2.)*cos((p*pi)/180.)**2.*&
              sin((p*pi)/180.)**2.)/&
              (L+2.*m)-(C11*C12*sin((p*pi)/180.)**4.)/(L+2.*m))*&
              sin((pi*t)/180.)**2.*tan((pi*t)/180.)**2.;
          
          
          
          F1222(t,p)=-((2.*C13*C66*cos((p*pi)/180.)*sin((p*pi)/180.)*&
              sin((pi*t)/180.)**2.)/&
              (L+2.*m))+cos((p*pi)/180.)*sin((p*pi)/180.)*&
              (-((2.*C12*C66*cos((p*pi)/180.)**2.)/(L+2.*m))-&
              (2.*C11*C66*sin((p*pi)/180.)**2.)/(L+2.*m))*sin((pi*t)/180.)**2.*&
              tan((pi*t)/180.)**2.;
          
          
          F2222(t,p)=-(C13**2./(4.*(L+2.*m)))+&
              (1./2.)*(-((C12*C13*cos((p*pi)/180.)**2.)/(L+2.*m))-&
              (C11*C13*sin((p*pi)/180.)**2.)/(L+2.*m)+C13**2./(2.*(L+2.*m)))*sin((pi*t)/180.)**2.+&
              (1./2.)*(-((C12**2.*cos((p*pi)/180.)**4.)/(2.*(L+2.*m)))-&
              (C11*C12*cos((p*pi)/180.)**2.*sin((p*pi)/180.)**2.)/(L+2.*m)-&
              (C11**2*sin((p*pi)/180.)**4.)/(2.*(L+2.*m)))*sin((pi*t)/180.)**2.*&
              tan((pi*t)/180.)**2.;
          
          Rpp_aniso(t,p)=.5*(Dexm)+.5*(((Ddxm-8.*(VsAve**2.)*Dgxm/(VpAve**2.)))*cos(p*pi/180.)**2.+(Ddym&
		-8.*(VsAve**2.)*Dgxm/(VpAve**2.))*sin(p*pi/180.)**2.+2.*(Dchizm-4.*(VsAve**2.)*De45m/(VpAve**2.))*&
		sin(p*pi/180.)*cos(p*pi/180.)-Dezm)*sin(t*pi/180.)**2. &
              +.5*(Dexm*cos(p*pi/180.)**4.+Deym*sin(p*pi/180.)**4.+Ddzm*(cos(p*pi/180.)**2.)*(sin(p*pi/180.)**2.)+&
		2.*(De16m*cos(p*pi/180.)**2.+De26m*sin(p*pi/180.)**2.)*sin(p*pi/180.)*cos(p*pi/180.))*sin(t*pi/180.)**2.*tan(t*pi/180.)**2.;
          
          Rpp_FRAC(t,p)=F11(t,p)*alpha11+F12(t,p)*alpha12+F22(t,p)*alpha22+F1111(t,p)*beta1111+&
		F1112(t,p)*beta1112+F1122(t,p)*beta1122+F1222(t,p)*beta1222+F2222(t,p)*beta2222;
          
          Rpp(t,p)=Rpp0(t)+Rpp_aniso(t,p)+Rpp_FRAC(t,p);
      end do
      
  end do
  
   ! write (0,*) 'dim F11', shape(F11)
    !write (0,*) ' F11', F11
!################################################################################################
  jjj=size(F11,2)
  hhh=size(F11,1)
  cc=hhh*jjj
  allocate (Rpp_remainder1(cc))
  allocate (Rpp_aniso1(cc))
  allocate (Rpp_FRAC1(cc))
  allocate (a1 (cc))
  allocate (b1 (cc))
  allocate (c1 (cc))
  allocate (d1 (cc))
  allocate (e1 (cc))
  allocate (f1 (cc))
  allocate (g1 (cc))
  allocate (h1 (cc))
  allocate (F1D (cc*8))
  allocate (FF (8,cc))




  a1=reshape(F11/m,(/cc/))
  b1=reshape(F12/m,(/cc/))
  c1=reshape(F22/m,(/cc/))
  d1=reshape(F1111/m,(/cc/))
  e1=reshape(F1112/m,(/cc/))
  f1=reshape(F1122/m,(/cc/))
  g1=reshape(F1222/m,(/cc/))
  h1=reshape(F2222/m,(/cc/))

  Rpp_aniso1=reshape(Rpp_aniso,(/cc/))
  Rpp_FRAC1=reshape(Rpp_FRAC,(/cc/))

  F1D=(/ a1,b1,c1,d1,e1,f1,g1,h1 /)    !matrix of F values

  FF=reshape(F1D,(/ 8,cc /), order=(/2,1/)) 			! reshape 2D to 3D




!%%%%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%    			       %%%%%%%%%%%%%%%%%%%%%%
  allocate (albet (8))

  albet=m*(/ alpha11,alpha12,alpha22,beta1111,beta1112,beta1122,beta1222,beta2222/)
            	 write (0,*) 'albet', albet

  Rpp0_background=spread(Rpp0,2,phi); !Rpp due to isotropic materials
  

	allocate (noise(cc))
	allocate (noise1(theta,phi))
	allocate (noise2(theta,phi))


	call randomv(cc,noise)
          	 write (0,*) 'MAx rand',  maxval(noise)
          	 write (0,*) 'Min rand',  minval(noise)
	noise2= reshape(noise,(/theta,phi/))

	noise1=noise2*abs((sum(Rpp)/(size(Rpp,1)*size(Rpp,2)))/SNR)
          	 write (0,*) 'Sum Rpp', abs((sum(Rpp)/(size(Rpp,1)*size(Rpp,2)))/SNR)
        !  	 write (0,*) 'noise1',  noise1(1:10,4:8)
	Rpp_noisy=Rpp+noise1;
          !	 write (0,*) 'noise1',  noise1
          	 write (0,*) 'Rpp(15,15)',  Rpp(15,15)
          	 write (0,*) 'FF(6,15)',  FF(6,15)
	Rpp_remainder=Rpp_noisy-Rpp0_background-Rpp_aniso; 
	Rpp_Minus_iso_background=Rpp-Rpp0_background;

	Rpp_Noisy_Minus_background=Rpp_noisy-Rpp0_background;

 	Rpp_remainder1=reshape(Rpp_remainder,(/cc/))
 !         	 write (0,*) 'noise2', Rpp_remainder1


!      write (0,*) 'dim FF', shape(FF)
!      write (0,*) 'dim F11', shape(F11)
!      write (0,*) 'dim Rpp_remainder1', shape(Rpp_remainder1)

 	 allocate (cons (8))
         cons=matmul(FF,Rpp_remainder1)

!        write (0,*) 'cons',  cons


	phi_s1_actual=(90./pi)*atan(2.*alpha12/(alpha11-alpha22))  ! actual fast shear direction


        allocate (FFT (8,8))
        FFT=matmul(FF,transpose(FF))	

	allocate (weight1(8))

        nu=8
        I=30
        call cg(nu,cons,FFT,tolerancer,itmax,I,weight1)
          	 write (0,*) 'weight1',  weight1
  !C_pred=C+weight
	allocate (w1d(16))
	allocate (w2d(8,2))

       W1d=(/ albet,weight1 /)    !matrix of F values

       W2d=reshape(W1d,(/ 8,2 /), order=(/2,1/)) 			! reshape 2D to 3D
          	 write (0,*) 'W2d',  shape(W2d)
!#################################################################################################################################


  !call rsf_write(commonoffset,rv20)   ! write out the common offset section
        call rsf_write(aw,albet)     ! write out the 3D volume
        call rsf_write(w1,weight1)     ! write out the 3D volume
        call rsf_write(refl,Rpp)     ! write out the 3D volume
        call rsf_write(NoisyRps,Rpp_noisy)     ! write out the 3D volume
        call rsf_write(Winvpred,Rpp_noisy)     ! write out the 3D volume
   
    
end program converted


