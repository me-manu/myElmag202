!=============================================================================!
!=============================================================================!
!                 inititialization of general purpose tools                   !
!=============================================================================!
!=============================================================================!
subroutine init(myid,n_proc)
!-----------------------------------------------------------------------------
! calls: banner,init_EBL,tabulate_rates,init_arrays
!-----------------------------------------------------------------------------
!f2py intent(in,out) myid
!f2py intent(in,out) n_proc
  use internal, only  : iseed
!  use user_result, only : filename
  use cosmology
  implicit none
  integer myid,n_proc

  if (myid==0) then
     call banner(n_proc,0)
     write(*,*) '!!!!!  init starts  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  end if

  iseed = 15321+2*myid      ! initialisation for randon number (NumRec)

  t_0 = 2.d0/(3.d0*sqrt(Omega_v)*H_0) * log((1.d0+sqrt(Omega_v))/sqrt(Omega_m))
  if (myid==0) then
     open(99,file='error')
!     write(*,*) 'cosmological parameters:'
!     write(*,*) 'Omega_m = ',Omega_m
!     write(*,*) 'Omega_v = ',Omega_v
!     write(*,*) 'age of the Universe/sec ',t_0
!     write(*,*) 
  end if

  call init_EBL(myid)           !initialize EBL
  call tabulate_rates           !tabulate integrated EBL density,
                                !m.f.p. of particles, ICS energy loss (E<Ethr)
  if (myid==0) write(*,*) '!  elmag. rates tabulated '
  call init_arrays(myid)        !initialize redshift tables

  if (myid==0)  then
     write(*,*) '!  init done '
!     write(*,*) '!  files will be saved *',filename
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) 
  end if

end subroutine init
!=============================================================================!
!                   init_EBL - initialization of EBL                          !
!=============================================================================!
!=============================================================================!
subroutine init_EBL(myid)
!-----------------------------------------------------------------------------
! calls: error
!-----------------------------------------------------------------------------
!f2py intent(in) myid
  use user_variables, only : model,tablefile_n,tablefile_z
  use EBL_fit  
  implicit none
  integer :: myid,i,j,k

  select case (model)
  case(1)                     ! 'best-fit' model of Kneiske & Doll
     n_EBL1 = 201
     allocate (z_ir(n_EBL1))
     allocate (eir(n_EBL1))
     allocate (anir(n_EBL1,n_EBL1))
     open(1,file=tablefile_n,status='old')
     do i=1,n_EBL1
        read (1,*) eir(i),anir(i,:)
     enddo
     close(1)     
       if (myid==0) write(*,*) '!  read in best fit' 

  case(2)                     ! 'lower limit' model of Kneiske & Doll
     n_EBL1 = 160
     n_EBL2 = 200
     allocate (z_ir(n_EBL2))
     allocate (eir(n_EBL1))
     allocate (anir(n_EBL1,n_EBL2))
     open(1,file=tablefile_z,status='old')
     do i=1,n_EBL2
        read(1,*)k,z_ir(i)
     enddo
     close(1)     
     open(1,file=tablefile_n,status='old')
     do i=1,n_EBL1
        read (1,*)eir(i),anir(i,:)
     enddo
     close(1)   
     if (myid==0) write(*,*) '!  read in lower limit' 
 
 case(3)                      ! Franceschini model, arXiv:0805.1841; data taken directly from paper
    n_EBL1 = 31
    n_EBL2 = 11
    allocate (z_ir(n_EBL2))
    allocate (eir(n_EBL1))
    allocate (anir(n_EBL1,n_EBL2)) 
    do i=1,n_EBL2
       z_ir(i)=0.2d0*(i-1)
    enddo    
    open(1,file=tablefile_n,status='old')
    do i=1,n_EBL1
       read(1,*)eir(i),anir(i,:)  
    enddo
    close(1)
    if (myid==0) write(*,*) '!  read in model from Franceschini et al. (2008) Franceschini' 
    
  case(4)                      ! Model C from Finke et al. (2010), arXiv:0905.1115; data taken from http://www.phy.ohiou.edu/~finke/EBL/
    n_EBL1 = 397
    n_EBL2 = 500
    allocate (z_ir(n_EBL2))
    allocate (eir(n_EBL1))
    allocate (anir(n_EBL1,n_EBL2)) 
    open(1,file=tablefile_n,status='old')
     do i=1,n_EBL2
       z_ir(i)=0.01d0*(i-1)
     enddo    
    do i=1,n_EBL1
       read(1,*)eir(i),anir(i,:)
    enddo
    close(1)
    if (myid==0) write(*,*) '!  read in model C from Finke et al. (2010)' 
  
  case(5)                     ! Model from Gilmore et al. (2012), arXiv:1104.0671; data taken from physics.ucsc.edu/~joel/EBLdata-Gilmore2012/
     n_EBL1 = 101
     n_EBL2 = 18
     allocate (z_ir(n_EBL2))
     allocate (eir(n_EBL1))
     allocate (anir(n_EBL1,n_EBL2))
     open(1,file=tablefile_z,status='old')
     do i=1,n_EBL2
        read(1,*)k,z_ir(i)
     enddo
     close(1)     
     open(1,file=tablefile_n,status='old')
     do i=1,n_EBL1
        read (1,*)eir(i),anir(i,:)
     enddo
     close(1)   
     if (myid==0) write(*,*) '!  read in model from Gilmore et al. (2012)' 

  case(6)                     ! Model from Dominguez et el.
     n_EBL1 = 50
     n_EBL2 = 17
     allocate (z_ir(n_EBL2))
     allocate (eir(n_EBL1))
     allocate (anir(n_EBL1,n_EBL2))

     open(1,file=tablefile_z,status='old')
     do i=1,n_EBL2
        read(1,*)k,z_ir(i)
     enddo
     close(1)     
     open(1,file=tablefile_n,status='old')
     do i=1,n_EBL1
        read (1,*)eir(i),anir(i,:)
     enddo
     close(1)   
     if (myid==0) write(*,*) '!  read in model from Dominguez et al. (2011)' 
  !End: New models
  
  case default
     call error('wrong EBL model',0)
  end select

end subroutine init_EBL
!=============================================================================!
! tabulation of integrated EBL density, m.f.p., ICS energy loss (E<Ethr)      !
!=============================================================================!
subroutine tabulate_rates 
!-----------------------------------------------------------------------------
! calls: w_EBL_density_tab,rate_EBL_tab,eloss_thr_tab
!-----------------------------------------------------------------------------
  use constants
  use EBL_fit
  use user_variables, only : ethr,egmax
  use xsec
  implicit none
  integer icq,k,l
  double precision e0,emin,zz   
  double precision rate_EBL_tab,w_EBL_density_tab,eloss_thr_tab

!    interaction with photon background (integrated background density)
  do k=1,50 
     if(k.le.11) then
        emin=1.d-15*(erad1*1.d15)**((k-1)/10.d0)  !E-cutoff for the background 
     else
        emin=erad1*(eirmax/erad1)**((k-11)/40.d0)
     end if
     do l=1,11
        zz=.5d0*(l-1)                             !redshift
        denc(k,l)=log(w_EBL_density_tab(emin,zz)*emin**(1.d0-brad1))
     end do
  end do
    
!    interaction with photon background (inverse m.f.p.)
  do icq=1,2
     if (icq.eq.1) then
        emin=ame**2/eirmax                        ! E-cutoff for HE photon
     else
        emin=ethr+.5d0*ame**2/eirmax &            ! E-cutoff for HE electron
       /(1.d0+dsqrt(1.d0+ame**2/ethr**2/eirmax*(ethr-eirmax)))
     end if
     do k=2,51
        e0=emin*(egmax/emin)**((k-1)/50.d0)       ! photon/electron energy
        do l=1,11
           zz=.5d0*(l-1)                          ! redshift
           if (k.eq.1) then       
              sigc(k,l,icq)=-80.d0
           else
              sigc(k,l,icq)=log(rate_EBL_tab(e0,zz,icq-1))   ! inverse m.f.p.
           end if	   	  
        end do
     end do
  end do
  
!    interaction with photon background (E-loss due to ICS under threshold)
  do k=2,51
     e0=ame*(egmax/ame)**((k-1)/50.d0)            ! electron energy
     do l=1,11
        zz=.5d0*(l-1)                             ! redshift
        zsigc(k,l)=log(eloss_thr_tab(e0,zz))      ! ICS energy-loss
     end do
  end do
 
end subroutine tabulate_rates
!=============================================================================!
!=============================================================================!
!       tabulation of interaction rate/cm on photon background                !
!=============================================================================!
double precision function rate_EBL_tab(e0,zz,icq)
!-----------------------------------------------------------------------------
! calls: w_EBL_density,sigics,sigpair
!-----------------------------------------------------------------------------
!f2py intent(in) e0
!f2py intent(in) zz
!f2py intent(in) icq
!f2py intent(out) rate_EBL_tab
  use constants, only : ame
  use EBL_fit
  use user_variables, only : ethr
  implicit none
  integer i,icq,m
  double precision e0,zz,sig1,sig2,sig3,sig4,de,emin,em,emin1,emax1,s,dsq,power
  double precision x1(7),a1(7)
  double precision w_EBL_density,sigics,sigpair
  data x1/.9862838d0,.9284349d0,.8272013d0,.6872929d0,.5152486d0,.3191124d0,.1080549d0/
  data a1/.03511946d0,.08015809d0,.1215186d0,.1572032d0,.1855384d0,.2051985d0,.2152639d0/
     
  rate_EBL_tab=0.d0
  sig1=0.d0
  sig2=0.d0
  sig3=0.d0
  sig4=0.d0

  select case (icq)
  case (0)                                                      ! photons

     de=4.d0*e0
     emin=ame**2/e0

     if (emin.ge.eirmax) then       !below the threshold for pair production
          rate_EBL_tab=0.d0
        return
     end if

     if (emin.lt.erad1) then
        emin1=emin**brad1
        emax1=erad1**brad1
        do i=1,7
           do m=1,2
             em=(.5d0*(emin1+emax1+(2*m-3)*x1(i)*(emax1-emin1)))**(1.d0/brad1)
             s=em*de
             sig1=sig1+a1(i)*em**(2.d0-brad1)*sigpair(s)*w_EBL_density(em,zz)
           enddo
        enddo
        sig1=sig1*(emax1-emin1)/brad1
     endif

     if (emin.lt.erad2)then
        do i=1,7
           do m=1,2
              em=2.d0/(1.d0/erad2+1.d0/max(emin,erad1) &
	      +(2*m-3)*x1(i)*(1.d0/max(emin,erad1)-1.d0/erad2))
              s=em*de
              sig2=sig2+a1(i)*em**3*sigpair(s)*w_EBL_density(em,zz)
           enddo
        enddo
        sig2=sig2*(1.d0/max(emin,erad1)-1.d0/erad2)
     endif

     if (emin.lt.eirmin*(1.d0+zz))then
        do i=1,7
           do m=1,2
              em=max(emin,erad2)*(eirmin*(1.d0+zz)/max(emin,erad2)) &
	      **(.5d0+x1(i)*(m-1.5d0))
              s=em*de
              sig3=sig3+a1(i)*em**2*sigpair(s)*w_EBL_density(em,zz)
           enddo
        enddo
        sig3=sig3*dlog(eirmin*(1.d0+zz)/max(emin,erad2))
     endif

     do i=1,7
        do m=1,2
           em=eirmax/sqrt(.5d0*((eirmax/max(emin,eirmin*(1.d0+zz)))**2 &
	   +1.d0+(2*m-3)*x1(i)*((eirmax/max(emin,eirmin*(1.d0+zz)))**2-1.d0)))
           s=em*de
           sig4=sig4+a1(i)*em**4*sigpair(s)*w_EBL_density(em,zz)
        enddo
     enddo
     sig4=sig4*(1.d0/max(emin,eirmin*(1.d0+zz))**2-1.d0/eirmax**2)/2.d0
     rate_EBL_tab=sig1+sig2+sig3+sig4
       
  case (-1,1)                                                      ! e^+(e^-)

     dsq=dsqrt(max(0.d0,1.d0-ame**2/e0**2))
     de=2.d0*e0*(1.d0+dsq)
     emin=ame**2*ethr/(e0-ethr)/de
     if (emin.ge.eirmax) then      !too low energy 
                                    !(both final particles below Ethr)
        rate_EBL_tab = 0.d0
        return             
     end if

     if (emin.lt.erad1)then
        emin1=emin**brad1
        emax1=erad1**brad1
        do i=1,7
           do m=1,2
             em=(.5d0*(emin1+emax1+(2*m-3)*x1(i)*(emax1-emin1))) &
	     **(1.d0/brad1)
             s=ame**2+em*de
	     sig1=sig1+a1(i)*sigics(e0,s)*w_EBL_density(em,zz)*em**(2.d0-brad1)
           enddo
        enddo
        sig1=sig1*(emax1-emin1)/brad1
     endif

     if (emin.lt.erad2)then
        if(erad2*de.gt.10.d0*ame**2)then
	 power=brad2
	else
	 power=brad2+1.d0
	endif
       emax1=erad2**power
        emin1=max(emin,erad1)**power
        do i=1,7
           do m=1,2
              em=(.5d0*(emin1+emax1+(2*m-3)*x1(i)*(emax1-emin1)))**(1.d0/power)
              s=ame**2+em*de
              sig2=sig2+a1(i)*sigics(e0,s)*w_EBL_density(em,zz)*em**(2.d0-power)
           enddo
        enddo
        sig2=sig2*(emax1-emin1)/power
     endif

     if (emin.lt.eirmin*(1.d0+zz))then
        emax1=exp(-eirmin/akt)
        emin1=exp(-max(emin,erad2)/akt/(1.d0+zz))
        do i=1,7
           do m=1,2
              em=-akt*(1.d0+zz)*dlog(.5d0*(emin1+emax1 &
	      +(2*m-3)*x1(i)*(emax1-emin1)))
              s=ame**2+em*de
              sig3=sig3+a1(i)*sigics(e0,s)*w_EBL_density(em,zz)*em &
	      *exp(em/akt/(1.d0+zz))
           enddo
        enddo
        sig3=-sig3*(emax1-emin1)*akt*(1.d0+zz)
     endif

     do i=1,7
        do m=1,2
           em=eirmax/dsqrt(.5d0*((eirmax/max(emin,eirmin*(1.d0+zz)))**2 &
	   +1.d0+(2*m-3)*x1(i)*((eirmax/max(emin,eirmin*(1.d0+zz)))**2-1.d0)))
           s=ame**2+em*de
           sig4=sig4+a1(i)*sigics(e0,s)*w_EBL_density(em,zz)*em**4
        enddo
     enddo
     sig4=sig4*(1.d0/max(emin,eirmin*(1.d0+zz))**2-1.d0/eirmax**2)/2.d0 
     rate_EBL_tab=(sig1+sig2+sig3+sig4)/4.d0/dsq*(1.d0+dsq)**2

  end select

end function rate_EBL_tab
!=============================================================================!
!=============================================================================!
!       tabulation of energy loss/cm for ICS (E<E_thr) on photon background   !
!=============================================================================!
double precision function eloss_thr_tab(e0,zz)
!-----------------------------------------------------------------------------
! calls: w_EBL_density,zsigics
!-----------------------------------------------------------------------------
!f2py intent(in) e0
!f2py intent(in) zz
!f2py intent(out) eloss_thr_tab
  use constants, only : ame
  use EBL_fit
  use user_variables, only : ethr
  implicit none
  integer i,m
  double precision e0,zz,sig1,sig2,sig3,sig4,de,dsq,emin,s,emin1,emax1
  double precision x1(7),a1(7)
  double precision w_EBL_density,zsigics
  data x1/.9862838d0,.9284349d0,.8272013d0,.6872929d0,.5152486d0,.3191124d0,.1080549d0/
  data a1/.03511946d0,.08015809d0,.1215186d0,.1572032d0,.1855384d0,.2051985d0,.2152639d0/
     
  eloss_thr_tab=0.d0
  sig1=0.d0
  sig2=0.d0
  sig3=0.d0
  sig4=0.d0

  dsq=dsqrt(max(0.d0,1.d0-ame**2/e0**2))
  de=2.d0*e0*(1.d0+dsq)

  do i=1,7
     do m=1,2
        emin=erad1*(.5d0+x1(i)*(m-1.5d0))**(1.d0/brad1)
        s=ame**2+emin*de
        sig1=sig1+a1(i)*zsigics(e0,s)*w_EBL_density(emin,zz)*emin**(2.d0-brad1)
     enddo
  enddo
  sig1=sig1*erad1**brad1/brad1

  emin1=erad1**brad2
  emax1=erad2**brad2
  do i=1,7
     do m=1,2
        emin=(.5d0*(emin1+emax1+(2*m-3)*x1(i)*(emax1-emin1)))**(1.d0/brad2)
        s=ame**2+emin*de
        sig2=sig2+a1(i)*zsigics(e0,s)*w_EBL_density(emin,zz)*emin**(2.d0-brad2)
     enddo
  enddo
  sig2=sig2*(emax1-emin1)/brad2

  do i=1,7
     do m=1,2
        emin=erad2*(eirmin*(1.d0+zz)/erad2)**(.5d0+x1(i)*(m-1.5d0))
        s=ame**2+emin*de
        sig3=sig3+a1(i)*zsigics(e0,s)*w_EBL_density(emin,zz)*emin**2
     enddo
  enddo
  sig3=sig3*dlog(eirmin*(1.d0+zz)/erad2)

  do i=1,7
     do m=1,2
        emin=eirmax/dsqrt(.5d0*((eirmax/eirmin/(1.d0+zz))**2+1.d0  &
	+(2*m-3)*x1(i)*((eirmax/eirmin/(1.d0+zz))**2-1.d0)))
        s=ame**2+emin*de
        sig4=sig4+a1(i)*zsigics(e0,s)*w_EBL_density(emin,zz)*emin**4
     enddo
  enddo
  sig4=sig4*(1.d0/(eirmin*(1.d0+zz))**2-1.d0/eirmax**2)/2.d0
  
  eloss_thr_tab=(sig1+sig2+sig3+sig4)/4.d0/dsq*(1.d0+dsq)**2

end function eloss_thr_tab  
!=============================================================================!
!=============================================================================!
!          tabulation of weighted background photon density Eq. (9)           ! 
!=============================================================================!
double precision function w_EBL_density_tab(emin,zz)
!-----------------------------------------------------------------------------
! calls: aintIR
!-----------------------------------------------------------------------------
!f2py intent(in) emin
!f2py intent(in) zz 
!f2py intent(out) w_EBL_density_tab
  use EBL_fit
  implicit none
  integer i,m
  double precision emin,zz,z,zmin,zmax,emax1,emin1,eg,aktz,denc,wbb_density
  double precision x1(7),a1(7)
  data x1/.9862838d0,.9284349d0,.8272013d0,.6872929d0,.5152486d0, &
 .3191124d0,.1080549d0/
  data a1/.03511946d0,.08015809d0,.1215186d0,.1572032d0, &
 .1855384d0,.2051985d0,.2152639d0/
  double precision aintIR

  wbb_density=0.d0
  w_EBL_density_tab = 0.d0
  if (emin.ge.eirmax) return              !threshold above IR-cutoff for EBL
  
  emin1=max(emin,eirmin*(1.d0+zz))
  do i=1,7
     do m=1,2
        eg=eirmax/(.5d0*((eirmax/emin1)**3+1.d0+(2*m-3)*x1(i) &
             *((eirmax/emin1)**3-1.d0)))**(1.d0/3.d0)
        wbb_density=wbb_density+a1(i)*eg**2*aintIR(eg,zz)
     end do
  end do

  wbb_density=wbb_density/6.d0*(1.d0/emin1**3-1.d0/eirmax**3)

  aktz=akt*(1.d0+zz)
  zmax=exp(-emin/aktz)
  wbb_density=wbb_density-aktz*1.32d13*dlog(1.d0-zmax)
      
  if (emin.lt.erad2) wbb_density=wbb_density+arad2/(brad2-1.d0)  &
  *(erad2**(brad2-1.d0)-max(emin,erad1)**(brad2-1.d0))
     
  if (emin.lt.erad1) then
     denc=0.d0
     if (emin.gt.0.d0) then
        emin1=emin**(brad1-1.d0)
     else
        emin1=0.d0
     end if
     emax1=erad1**(brad1-1.d0)
     do i=1,7
        do m=1,2
           eg=(.5d0*(emax1+emin1+(2*m-3)*x1(i)*(emax1-emin1)))**(1.d0/(brad1-1.d0))
           denc=denc+a1(i)/(1.d0+crad*eg**drad)
        end do
     end do
     wbb_density=wbb_density+denc*arad1/(brad1-1.d0)/2.d0*(emax1-emin1)
  end if

  w_EBL_density_tab = wbb_density

end function w_EBL_density_tab
!=============================================================================!
!                    EBL density (interpolation)                              ! 
!=============================================================================!
double precision function aintIR(E,z)
!f2py intent(in) E
!f2py intent(in) z
!f2py intent(out) aintIR
  use user_variables, only : model,ir_rescale
  use EBL_fit
  implicit none
  integer ie,iz,je,jz,jzmax,jeF
  double precision e,z,zz,we(2),wz(2),eminF,emaxF,rparam,tparam

  aintIR=0.d0

  select case (model)
  case (1)
     if (e.ge.eir(201).or.e.le.eir(1)) return          
     
     if (z.le.1.d-3) then
        jz=1
        wz(1)=1.d0
        wz(2)=0.d0
        jzmax=1
     elseif (z.le..025d0) then
        jz=1
        wz(2)=(z-1.d-3)/.024d0
        wz(1)=1.d0-wz(2)
        jzmax=2
     else
        jzmax=2
        zz=z/0.025d0+1.d0
        jz=min(int(zz),200)
        wz(2)=zz-jz
        wz(1)=1.d0-wz(2)
     endif
     
     do je=1,200                                                 ! n_EBL-1
        if (e.gt.eir(je).and.e.le.eir(je+1)) goto 1
     enddo
     stop'wrong je'
1    we(2)=dlog(e/eir(je))/log(eir(je+1)/eir(je)) 
     we(1)=1.d0-we(2)
      
     do ie=1,2
        do iz=1,jzmax
           aintIR=aintIR+log(anir(je+ie-1,jz+iz-1))*wz(iz)*we(ie)
        enddo
     enddo

  case (2)
     if (e.lt.eir(160).or.e.ge.eir(1)) return
     
     if(z.le.z_ir(1))then
        jz=1
        wz(1)=1.d0
        wz(2)=0.d0
        jzmax=1
     else
        jzmax=2
        do jz=1,199
           if(z.gt.z_ir(jz).and.z.le.z_ir(jz+1))goto 11
        enddo
        stop'wrong jz'
11      jz=min(jz,198)
        wz(2)=(z-z_ir(jz))/(z_ir(jz+1)-z_ir(jz))
        wz(1)=1.d0-wz(2)
     endif
      
     do je=1,159
        if(e.lt.eir(je).and.e.ge.eir(je+1))goto 22
     enddo
     stop'wrong je'
22   we(2)=dlog(e/eir(je))/log(eir(je+1)/eir(je)) 
     we(1)=1.d0-we(2)
     
     do ie=1,2
        do iz=1,jzmax
           aintIR=aintIR+log(anir(je+ie-1,jz+iz-1))*wz(iz)*we(ie)
      enddo
      enddo
  
case (3)
     !if (e.lt.eir(31).or.e.ge.eir(1)) return
     if (e.lt.eir(31).or.e.ge.eir(1)) then
       aintIR=-15.d0
       goto 33
     endif     
     
     if(z.le.z_ir(1))then
        jz=1
        wz(1)=1.d0
        wz(2)=0.d0
        jzmax=1
     else
        jzmax=2
        do jz=1,10
           if(z.gt.z_ir(jz).and.z.le.z_ir(jz+1))goto 31
        enddo
        if(z.gt.2.d0) then
          jz=10
          goto 31
        endif
        stop'wrong jz'
31      jz=min(jz,9)
        wz(2)=(z-z_ir(jz))/(z_ir(jz+1)-z_ir(jz))
        wz(1)=1.d0-wz(2)
     endif
     
     do je=1,30
        if(e.lt.eir(je).and.e.ge.eir(je+1))goto 32
     enddo
     stop'wrong je'
32   we(2)=dlog(e/eir(je))/log(eir(je+1)/eir(je)) 
     we(1)=1.d0-we(2)

     do ie=1,2
        do iz=1,jzmax
           aintIR=aintIR+log(anir(je+ie-1,jz+iz-1))*wz(iz)*we(ie)
        enddo
     enddo
33   aintIR=aintIR
    
case (4)
     if (e.lt.eir(397).or.e.ge.eir(1)) return
     
     if(z.le.z_ir(1))then
        jz=1
        wz(1)=1.d0
        wz(2)=0.d0
        jzmax=1
     else
        jzmax=2
        do jz=1,499
           if(z.gt.z_ir(jz).and.z.le.z_ir(jz+1))goto 41
        enddo
        if(z.eq.5.d0) then
          jz=499
          goto 41
        endif
        stop'wrong jz'
41      jz=min(jz,498)
        wz(2)=(z-z_ir(jz))/(z_ir(jz+1)-z_ir(jz))
        wz(1)=1.d0-wz(2)
     endif
      
     do je=1,396
        if(e.lt.eir(je).and.e.ge.eir(je+1))goto 42
     enddo
     stop'wrong je'
42   we(2)=dlog(e/eir(je))/log(eir(je+1)/eir(je)) 
     we(1)=1.d0-we(2)
     
     do ie=1,2
        do iz=1,jzmax
           aintIR=aintIR+log(anir(je+ie-1,jz+iz-1))*wz(iz)*we(ie)
     enddo
    enddo

case (5)
     if (e.lt.eir(101).or.e.ge.eir(1)) return
     
     if(z.le.z_ir(1))then
        jz=1
        wz(1)=1.d0
        wz(2)=0.d0
        jzmax=1
     else
        jzmax=2
        do jz=1,17
           if(z.gt.z_ir(jz).and.z.le.z_ir(jz+1))goto 51
        enddo
        stop'wrong jz'
51      jz=min(jz,16)
        wz(2)=(z-z_ir(jz))/(z_ir(jz+1)-z_ir(jz))
        wz(1)=1.d0-wz(2)
     endif
      
     do je=1,100
        if(e.lt.eir(je).and.e.ge.eir(je+1))goto 52
     enddo
     stop'wrong je'
52   we(2)=dlog(e/eir(je))/log(eir(je+1)/eir(je)) 
     we(1)=1.d0-we(2)
     
     do ie=1,2
        do iz=1,jzmax
           aintIR=aintIR+log(anir(je+ie-1,jz+iz-1))*wz(iz)*we(ie)
      enddo
      enddo

case (6)
     if (e.lt.eir(50).or.e.ge.eir(1)) return
     
     if(z.le.z_ir(1))then
        jz=1
        wz(1)=1.d0
        wz(2)=0.d0
        jzmax=1
     else
        jzmax=2
        do jz=1,16
           if(z.gt.z_ir(jz).and.z.le.z_ir(jz+1))goto 61
        enddo
        if(z.gt.3.9d0) then
          jz=16
          goto 61
        endif
        stop'wrong jz'
61      jz=min(jz,15)
        wz(2)=(z-z_ir(jz))/(z_ir(jz+1)-z_ir(jz))
        wz(1)=1.d0-wz(2)
     endif
      
     do je=1,49
        if(e.lt.eir(je).and.e.ge.eir(je+1))goto 62
     enddo
     stop'wrong je'
62   we(2)=dlog(e/eir(je))/log(eir(je+1)/eir(je)) 
     we(1)=1.d0-we(2)
     
     do ie=1,2
        do iz=1,jzmax
           aintIR=aintIR+log(anir(je+ie-1,jz+iz-1))*wz(iz)*we(ie)
      enddo
      enddo
  

  case default
     stop 'wrong EBL model'

  end select

  aintIR=exp(aintIR)*(1.d0+z)**3*ir_rescale
  
end function aintIR
!=============================================================================!
! initializion of redshift tables                                             !
!=============================================================================!
subroutine init_arrays(myid)
!f2py intent(in) myid
  use user_variables, only : tablefile_redshift
  implicit none
  integer i,n,n1,myid
  parameter (n=310,n1=69)
  double precision t1(0:n),z1(0:n),t2(0:n),z2(0:n),t3(0:n),r3(0:n), &
 &  z4(0:n),r4(0:n) 
  common /data_div1/ t1,z1
  common /data_div2/ t2,z2
  common /data_div3/ t3,r3
  common /data_div4/ z4,r4
 
  open(1,file=tablefile_redshift, status = 'old')
  do i=0,n
     read(1,*) t1(i),z1(i),r3(i)
     z2(n-i) = z1(i)         
     t2(n-i) = t1(i)
  enddo
  close(1)
  t3 = t1
  z4 = z1
  r4 = r3
  if (myid==0) write(*,*) '!  redshift table read'

end subroutine init_arrays
!=============================================================================!
!=============================================================================!
subroutine banner(n_proc,i)
!f2py intent(in) n_proc
!f2py intent(in) i
  use user_result
  implicit none
  integer n_proc,i
  character(8)  :: date
  character(10) :: time

  write(*,*) 
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  if (i==0) then
     write(*,*) '!  ELMAG v2.02                                           !'
     write(*,*) '!  Authors: M.Kachelriess, S. Ostapchenko, R.Tomas       !'
     write(*,*) '!  Refs.:  Comput.Phys.Commun. 183 (2012) 1036-1043      !'
     write(*,*) '!          arxiv[1106.5508]                              !'
     write(*,*) '!  see http://elmag.sourceforge.net for more info        !'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) '!  EBL models 3-5 contributed by Andrey Saveliev         !'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) '! Disclaimer: this program comes without any guarantees! !' 
     write(*,*) '! Beware of errors and use common sense                  !' 
     write(*,*) '! interpreting results!                                  !'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     call date_and_time(DATE=date,TIME=time)
     print *, '! The run starts on ',  date(7:8),'/', date(5:6),'/',date(1:4), &
            ' at ',time(1:2), ':', time(3:4),'                  !'
     write(*,*) '!  on',n_proc,' processor(s)                          !'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) '!  ELMAG  -- start of the main program                   !'
  else
     write(*,*) 
!     write(*,*) ' files saved as *',filename
     write(*,*) 
     call date_and_time(DATE=date,TIME=time)
     print *, '! The run finished on ',  date(7:8),'/', date(5:6),'/',date(1:4), &
            ' at ',time(1:2), ':', time(3:4)
     write(*,*) 
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) '!  Elmag - end of the main program                       !'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  end if
               
end subroutine banner
!=============================================================================!
!=============================================================================!
