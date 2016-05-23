!=============================================================================!
!=============================================================================!
!          subroutines handling cascade evolution and stack                   !
!=============================================================================!
!=============================================================================!
! (1+1)-dimensional treatment of e/m cascades in extragalactic space          !
!=============================================================================!
subroutine cascade(icq,e00,weight0,z_in)
!-----------------------------------------------------------------------------
! calls: int_length, interaction, eloss, zloss, get_particle, angle_delay,  
!        registerresult, z_t, t_z, bemf, themf, error
! input:
!        icq - initial particle type (0 - photon, +/-1 - electron/positron);
!        e00 - initial particle energy (eV);
!        weight0 - initial particle weight;
!        z_in - initial particle redshift position
!-----------------------------------------------------------------------------
!f2py intent(in) icq
!f2py intent(in) e00
!f2py intent(in) weight0 
!f2py intent(in) z_in 
  use constants
  use stack, only: jcmb
  use cosmology
  use user_variables, only : ethr,cohlnth
  use internal, only : rcmb,bal
  implicit none
  integer icq,ncohl,ierr,noint
  double precision e00,z_in,weight0
  double precision e0,x,zz,t,weight,the1,the2,xxc,xx,dt,begmf,theta,x0,z1,zz0
  double precision beta_red,t_in,de,dthet,dethr,zloss
  double precision t_z,z_t,int_length,themf,eloss,bemf

! initialization
  e0=e00
  weight=weight0

  x0=0.d0 
  zz0=z_in
  t = t_z(z_in)
  t_in = t*t_0
  rcmb = (1.d0-t)*t_0*3.d10             !light-travel time from the source (cm)

  theta=0.d0                      
  the1=0.d0 
  the2=0.d0 
  xxc=0.d0
  xx=0.d0
  dt=0.d0                   

  bal=bal+e0*weight
  jcmb=0
  noint=0

  do
  
! new particle from stack
     if (jcmb.ne.0.and.noint.eq.0) call get_particle(e0,x0,zz0,t,weight,the1,the2,xxc,xx,dt,icq)
     if (.not.e0.gt.0.d0) call error('strange e0 in cascade',0)
     
     noint=0
  
     x=min(rcmb,int_length(e0,x0,zz0,icq))      ! interaction point
     if (x<x0) call error('x<x0 in cascade',0)
      
     t = t+(x-x0)/(t_0*3.d10)

     zz = z_t(t)                                ! new redshift
     
     if(zz0-zz.gt..05d0)then                    ! dz_max = 0.05 
      zz=zz0-.05d0
      noint=1
      t=t_z(zz)
      x=x0+(t-t_z(zz0))*t_0*3.d10
     endif
     
     if (zz.lt.0.d0) then
        if (zz.lt.-1.d-8) call error('zz<0 in cascade',0)
        zz=0.d0
     end if
     z1 = 1.d0+zz
     beta_red = e0 * H_0 * sqrt( Omega_m*z1**3.d0 + Omega_v )
     e0 = e0-(x-x0)/3.d10 * beta_red
      
     if (icq.ne.0) then                          ! electron/positron
        begmf= bemf(rcmb-x) 
        de=min(e0,eloss(e0,begmf)*(x-x0))        ! local magnetic field
        dethr=min(e0-de,e0*zloss(e0,zz0)*(x-x0)) ! synchrotron E-loss
        dthet=themf(e0,begmf)                    ! deflection angle per cm
        if (x-x0+xxc.lt.cohlnth) then
           the1=the1+dthet*(x-x0)                ! accumulate deflection angle
           xxc=xxc+x-x0
        else
           the1=the1+dthet*(cohlnth-xxc)         ! accumulate deflection angle
           the2=the2+the1**2                     ! accumulate deflection angle^2
           xxc=x-x0+xxc-cohlnth
           if (xxc.ge.cohlnth) then              ! account for random B-field        
              ncohl=int(xxc/cohlnth)
              the2=the2+(dthet*cohlnth)**2*ncohl ! accumulate deflection angle^2
              xxc=xxc-cohlnth*ncohl
           end if
           the1=dthet*xxc                        ! accumulate deflection angle
        end if
        dt=dt+(x-x0)/cy*(ame/e0)**2/2.d0         ! kinematic time delay 
        xx=x                                     ! distance from the source
        bal=bal-(de+dethr)*weight
        e0=e0-de-dethr                           ! E of the current particle
     end if
     
     if (e0<=ethr.or.x.ge.0.999d0*rcmb) then     ! end of particle tracking
        if (icq.eq.0) call angle_delay(the2,xx,rcmb,theta,dt,e0) ! obs. angle
        if (e0>ethr) then
           call registerresult(e0,theta,dt,weight,icq) ! record final particle 
        else
           bal=bal-e0*weight
        end if
	noint=0
     elseif(noint.eq.0)then                         ! interaction with EBL
        call interaction(e0,x,zz,t,weight,the1,the2,xxc,xx,dt,icq,ierr)
        if (ierr==1) bal=bal-e0*weight
     else
        x0=x
	zz0=zz
     end if
     
     if (noint.eq.0.and.jcmb==0) exit                           ! empty stack
      
  end do
      
end subroutine cascade
!=============================================================================!
!=============================================================================!
! calculation of deflection angles/ time delays for observed gammas           !
!=============================================================================!
subroutine angle_delay(the2,xx,rcmb,theta,dt,e0)
!-----------------------------------------------------------------------------
! calls: psran, error
! input:
!        the2 - squared deflection angle for the cascade branch;
!        xx - distance from the source (cm) for the parent electron/positron;
!        rcmb - distance between the source and the observer;
!        e0 - photon energy
! output:
!        theta - observation angle (degree) for the photon;
!        dt - time delay (yr) for the photon
!-----------------------------------------------------------------------------
!f2py intent(in) the2
!f2py intent(in) xx
!f2py intent(in) rcmb
!f2py intent(in,out) e0
!f2py intent(out) theta
!f2py intent(out) dt
  use user_variables, only : th_jet
  use cosmology, only : t_0
  use constants, only : cy,pi
  implicit none
  double precision the2,theta,xx,rcmb,dt,e0,sbeta,cbeta,beta,psran

  if (the2>pi**2/4.d0) then     !beyond small-angle approx. -> take isotropic
     dt=t_0
     cbeta=2.d0*psran()-1.d0 
     beta=acos(cbeta)*180.d0/pi
     sbeta=sqrt(1.d0-cbeta**2)
     theta=asin(xx/rcmb*sbeta)*180.d0/pi                  ! observation angle
     if (beta-theta.gt.th_jet) e0=0.d0                    ! dismiss photon
     return
  elseif (the2.eq.0.d0) then
     theta=0.d0
  else
     sbeta=sin(sqrt(the2))
     theta=asin(xx/rcmb*sbeta)*180.d0/pi                  ! observation angle
     if (sqrt(the2)*180.d0/pi-theta.gt.th_jet) e0=0.d0    ! dismiss photon
     dt=dt+2.d0*xx/cy*(1.d0-xx/rcmb)*sbeta**2             ! time delay
  end if

  if (theta<0.d0) call error('theta<0 in angle_delay',1)

end subroutine angle_delay
!=============================================================================!
!=============================================================================!
! interaction with background photons                                         !
!=============================================================================!
subroutine interaction(e0,x0,zz,t,weight,the1,the2,xxc,xx,dt,icq,ierr)
!-----------------------------------------------------------------------------
! calls: sample_photon, sample_electron, zpair, zics, store_particle, error
! input:
!        e0 - particle energy;
!        x0 - particle distance from the source;
!        zz - particle redshift position;
!        weight - particle weight;
!        t - time;
!        the1 - accumulated deflection angle within B-field coh. length (for e+-);
!        the2 - squared deflection angle for the particle;
!        xxc - travel distance within B-field coherence length (for e+-);
!        xx - distance from the source (cm) for last e+- in the cascade branch;
!        dt - particle time delay;
!        icq - particle type (0 - photon, +/-1 - electron/positron)
! output:
!        ierr - error code (0 - o.k., 1 - error)
!-----------------------------------------------------------------------------
!f2py intent(in) e0
!f2py intent(in) x0
!f2py intent(in) zz
!f2py intent(in) weight
!f2py intent(in) t
!f2py intent(in) the1
!f2py intent(in) the2
!f2py intent(in) xxc
!f2py intent(in) xx
!f2py intent(in) dt
!f2py intent(in) icq
!f2py intent(out) ierr
  implicit none
  integer icq,ierr
  double precision e0,x0,zz,t,weight,the1,the2,xxc,xx,dt,z,sgam
  double precision zpair,zics
  
  ierr=0   
  select case (icq) 
  case (0)                                   ! pair production on EBL
     call sample_photon(e0,zz,sgam,ierr)     ! sample c.m. energy for interaction
     if (ierr==1) return
     z=zpair(sgam)                           ! E-share taken by electron
     call store_particle(z*e0,x0,zz,t,z,weight,0.d0,the2,0.d0,xx,dt,1) ! record e- 
     call store_particle((1.d0-z)*e0,x0,zz,t,(1.d0-z),weight  &        ! record e+
 &   ,0.d0,the2,0.d0,xx,dt,-1)
  case (-1,1)                                ! ICS on EBL
     call sample_electron(e0,zz,sgam,ierr)     
     if (ierr==1) return
     z=zics(e0,sgam)                         ! E-share taken by electron/positron
     call store_particle(z*e0,x0,zz,t,z,weight,the1,the2,xxc,xx,dt,icq) ! e+-
     call store_particle((1.d0-z)*e0,x0,zz,t,(1.d0-z),weight  &         !gamma
 &   ,0.d0,the2+the1**2,0.d0,xx,dt,0)
  case default
     call error('wrong icq in interaction',0)
  end select

end subroutine interaction
!=============================================================================!
!=============================================================================!
!        add a particle to stack                                              !
!=============================================================================!
subroutine store_particle(e0,x0,zz,t,ze,weight,the1,the2,xxc,xx,dt,icq)
!-----------------------------------------------------------------------------
! calls: psran, error
! input:
!        e0 - particle energy;
!        x0 - particle distance from the source;
!        zz - particle redshift position;
!        weight - particle weight;
!        t - time;
!        ze - share of the parent particle energy;
!        the1 - accumulated deflection angle within B-field coh. length (for e+-);
!        the2 - squared deflection angle for the particle;
!        xxc - travel distance within B-field coherence length (for e+-);
!        xx - distance from the source (cm) for last e+- in the cascade branch;
!        dt - particle time delay;
!        icq - particle type (0 - photon, +/-1 - electron/positron)
!-----------------------------------------------------------------------------
!f2py intent(in) e0
!f2py intent(in) x0
!f2py intent(in) zz
!f2py intent(in) weight
!f2py intent(in) t
!f2py intent(in) ze
!f2py intent(in) the1
!f2py intent(in) the2
!f2py intent(in) xxc
!f2py intent(in) xx
!f2py intent(in) dt
!f2py intent(in) icq
  use stack
  use internal, only : bal
  use user_variables, only : ethr,a_smp
  implicit none
  integer icq
  double precision e0,x0,zz,t,ze,weight,the2,the1,xxc,xx,dt
  integer i,n
  double precision aweight,psran

  if (e0.le.ethr) then                ! disregard underthreshold particles
     bal=bal-e0*weight
     return
  end if
      
  aweight=ze**a_smp                   ! sampling weight for produced particle
  if (psran().gt.aweight) return      ! disregard particle with prob. (1-aweight)

  jcmb=jcmb+1                         ! enlarge stack
  if (jcmb==n_max) call error('enlarge storage!',0)

  act = one_event(icq,e0,x0,zz,t,weight,the1,the2,xxc,xx,dt)

  if (jcmb.gt.1) then                 ! add particle to E-ordered stack
     do i=1,jcmb-1
        n=jcmb-i
        if(event(n)%en.gt.e0)then
           event(n+1)=act                     ! add particle to stack
           event(n+1)%w=event(n+1)%w/aweight  ! total weight of the particle
           return
        else
           event(n+1)=event(n)                ! re-arrange stack
        endif
     enddo
  endif
  event(1)=act                        ! 1st particle in stack
  event(1)%w=event(1)%w/aweight       ! total weight of the particle

end subroutine store_particle
!=============================================================================!
!=============================================================================!
!        get a particle from stack                                            !
!=============================================================================!
subroutine get_particle(e0,x0,zz,t,weight,the1,the2,xxc,xx,dt,icq)
!-----------------------------------------------------------------------------
! output:
!        e0 - particle energy;
!        x0 - particle distance from the source;
!        zz - particle redshift position;
!        weight - particle weight;
!        t - time;
!        the1 - accumulated deflection angle within B-field coh. length (for e+-);
!        the2 - squared deflection angle for the particle;
!        xxc - travel distance within B-field coherence length (for e+-);
!        xx - distance from the source (cm) for last e+- in the cascade branch;
!        dt - particle time delay;
!        icq - particle type (0 - photon, +/-1 - electron/positron)
!-----------------------------------------------------------------------------
!f2py intent(out) e0
!f2py intent(out) x0
!f2py intent(out) zz
!f2py intent(out) weight
!f2py intent(out) t
!f2py intent(out) the1
!f2py intent(out) the2
!f2py intent(out) xxc
!f2py intent(out) xx
!f2py intent(out) dt
!f2py intent(out) icq
  use stack
  implicit none
  integer icq
  double precision e0,x0,zz,t,weight,the1,the2,xxc,xx,dt

  icq=event(jcmb)%icq
  e0=event(jcmb)%en
  x0=event(jcmb)%x
  zz=event(jcmb)%z
  t=event(jcmb)%t
  weight=event(jcmb)%w
  the1=event(jcmb)%the1
  the2=event(jcmb)%the2
  xxc=event(jcmb)%xxc
  xx=event(jcmb)%xx
  dt=event(jcmb)%dt

  jcmb=jcmb-1

end subroutine get_particle
!=============================================================================!
!=============================================================================!
!            subroutines handling single interaction                          !
!=============================================================================!
!=============================================================================!
! sample c.m. energy for gamma-gamma_EBL interaction                          !
!=============================================================================!
subroutine sample_photon(e0,zz,sgam,ierr)
!-----------------------------------------------------------------------------
! calls: w_EBL_density, sigpair, psran
! input:
!        e0 - photon energy;
!        zz - photon redshift position
! output:
!        sgam - c.m. energy for gamma-gamma interaction;
!        ierr - error code (0 - o.k., 1 - error)
!
! NB: rejection method used - optimized for the EBL models of Kneiske & Doll;
! the procedure may have to be re-adjusted for different (higher EBL) models
!-----------------------------------------------------------------------------
!f2py intent(in) e0
!f2py intent(in) zz
!f2py intent(out) sgam
!f2py intent(out) ierr
  use EBL_fit
  use constants
  use user_variables, only : ethr,model
  use internal, only : debug
  implicit none
  integer ierr,nrej,nrejmax,nrenorm
  double precision e0,zz,sgam,de,emin,emin0 &
  ,etrans1,etrans2,aw1,aw2,aw3,gb,gb0,gbmax,gbnorm,rrr,gnorm1,gnorm2,gnorm3
  double precision psran,sigpair,w_EBL_density

  de=4.d0*e0  
  emin0=ame**2/e0                  ! minimal required energy for EBL photon
  if (emin0.ge.eirmax) then        ! wrong kinematics
!!!     write(*,*)'photon:',emin0,eirmax,e0
     ierr=1
     return
  end if
  
  nrej=0
  nrejmax=3000                     ! user-defined limit on the N of rejections
  nrenorm=0
  gbmax=0.d0
  gbnorm=2.5d0                      ! normalization factor for rejection

  etrans1=1.d-6                    ! parameters for 'proposal function'
  etrans2=eirmin*(1.d0+zz)**1.25

! partial weights for different energy intervals for 'proposal function'  
  if(emin0.lt.etrans1)then
   gnorm1=w_EBL_density(etrans1,zz)*etrans1*sigpair(etrans1*de)
   aw1=gnorm1*(etrans1/emin0-1.d0)*etrans1
  else
   aw1=0.d0
  endif

  if(emin0.lt.etrans2)then
   gnorm2=w_EBL_density(max(etrans1,2.d0*emin0),zz)*max(etrans1,2.d0*emin0) &
   *sigpair(max(etrans1,2.d0*emin0)*de)
   aw2=gnorm2*(exp((max(etrans1,2.d0*emin0)-max(etrans1,emin0)) &
   /akt/(1.d0+zz)/1.5)-exp((max(etrans1,2.d0*emin0)-etrans2) &
   /akt/(1.d0+zz)/1.5))*akt*(1.d0+zz)*1.5d0
  else
   aw2=0.d0
  endif
  
  gnorm3=w_EBL_density(eirmax/2.d0,zz)*eirmax*sigpair(eirmax*de)
  aw3=gnorm3*((eirmax/max(emin0,etrans2))**2.5d0-1.d0)/2.5d0*eirmax

1 rrr=psran()*(aw1+aw2+aw3)

! sample emin (= sgam/de) according to the 'proposal function';
! define 'rejection function' ( gb = f(emin) / f_proposal(emin) )
  if(rrr.lt.aw1)then
   emin=etrans1/(1.d0+psran()*(etrans1/emin0-1.d0))
   gb0=gnorm1*(etrans1/emin)**2
   
  elseif(rrr.lt.aw1+aw2)then
   emin=etrans2-akt*(1.d0+zz)*1.5d0*dlog(1.d0-psran() &
   *(1.d0-exp((etrans2-max(emin0,etrans1))/akt/(1.d0+zz)/1.5d0)))
   gb0=gnorm2*exp((max(etrans1,2.d0*emin0)-emin)/akt/(1.d0+zz)/1.5d0)
   
  else
   emin=eirmax/(1.d0-psran()*(1.d0-(eirmax/max(emin0,etrans2))**2.5d0))**.4d0
   gb0=gnorm3*(eirmax/emin)**3.5d0
  endif
  gb0=gb0*gbnorm

  if(model.eq.3)then
   if(emin0/etrans2.gt.1.d0)gb0=gb0*1.d2
   if(emin0/etrans2.gt.1.d2)gb0=gb0*3.d0
  endif
    
  sgam=emin*de                 ! c.m. energy for gamma-gamma interaction
  gb=w_EBL_density(emin,zz)*sigpair(sgam)*emin/gb0

!  if (gb.gt.1.d0.and.nrenorm.eq.0) write(*,*)'sample_cmb(photon): gb=' &
!  ,gb,nrenorm,emin,emin0,emin0/etrans2/1.d3
    
  if (psran().gt.gb) then       ! rejection
     nrej=nrej+1                ! total number of rejections for current sampling
     gbmax=max(gbmax,gb)        ! maximal value for rejection function
     if(nrej.gt.nrejmax)then    ! too many rejections
      if(gbmax.le.0.d0)then     ! wrong kinematics
       write(*,*)'photon: gbmax=0!!!'
       ierr=1
       return
      else
!       write(*,*)'nrej(gamma)>nrejmax',nrej,emin0/etrans2,nrenorm,e0/1.d12,gbmax
       gbnorm=gbnorm*gbmax*2.d0 ! change normalization for the rejection function
       gbmax=0.d0
       nrenorm=nrenorm+1
       nrej=0
      endif
     endif
     goto 1                     ! new try
  end if
  
end subroutine sample_photon

!=============================================================================!
!=============================================================================!
! sample c.m. energy for e(+-)-gamma_EBL interaction                          !
!=============================================================================!
subroutine sample_electron(e0,zz,sgam,ierr)
!-----------------------------------------------------------------------------
! calls: w_EBL_density, sigics, psran
! input:
!        e0 - electron/positron energy;
!        zz - electron/positron redshift position
! output:
!        sgam - c.m. energy for e(+-)-gamma interaction;
!        ierr - error code (0 - o.k., 1 - error)
!
! NB: rejection method used - optimized for the EBL models of Kneiske & Doll;
! the procedure may have to be re-adjusted for different (higher EBL) models
!-----------------------------------------------------------------------------
!f2py intent(in) e0
!f2py intent(in) zz
!f2py intent(out) sgam
!f2py intent(out) ierr
  use EBL_fit
  use constants
  use user_variables, only : ethr,model
  use internal, only : debug
  implicit none
  integer nrej,ierr,nrejmax,nrenorm
  double precision e0,zz,sgam,de,emin,emin0,gb,gb0,gbmax,gbnorm   &
& ,etrans1,etrans2,etrans3,aw1,aw2,aw3,aw4,rrr,gnorm1,gnorm2,gnorm3,gnorm4
  double precision psran,sigics,w_EBL_density

  de=2.d0*e0*(1.d0+dsqrt(max(0.d0,1.d0-ame**2/e0**2)))
  emin0=ame**2*ethr/(e0-ethr)/de    ! minimal required energy for EBL photon
  
  if (emin0.ge.eirmax) then         ! wrong kinematics
!!!     write(*,*)'electron:',emin0,eirmax,e0
     ierr=1
     return
  end if
  
  nrej=0
  nrejmax=3000                      ! user-defined limit on the N of rejections
  nrenorm=0
  gbmax=0.d0
  gbnorm=3.d0                       ! rejection normalization factor
  
  etrans1=2.d-10                    ! parameters for 'proposal function'
  etrans2=4.d-7
  etrans3=eirmin*(1.d0+zz)**1.2
  
! partial weights for different energy intervals for 'proposal function'  
  if(emin0.lt.etrans1)then
   gnorm1=emin0*w_EBL_density(emin0,zz)*sigics(e0,ame**2+2.d0*emin0*de)
   aw1=gnorm1*((etrans1/emin0)**(brad1+1.d0)-1.d0)/(brad1+1.d0)*emin0
  else
   aw1=0.d0
  endif

  if(emin0.lt.etrans2)then
   gnorm2=etrans2*w_EBL_density(etrans2,zz)*sigics(e0,ame**2+etrans2*de)
   if(etrans2*de.gt.10.d0*ame**2)then
    gnorm2=max(gnorm2,w_EBL_density(etrans1,zz)*sigics(e0,ame**2 &
    +etrans1*de)*etrans1*(etrans1/etrans2)**(1.d0-brad2))
    aw2=gnorm2*(1.d0-(max(emin0,etrans1)/etrans2)**brad2)/brad2*etrans2
   else
    gnorm2=max(gnorm2,w_EBL_density(etrans1,zz) &
    *sigics(e0,ame**2+etrans1*de)*etrans1*(etrans2/etrans1)**brad2)
    aw2=gnorm2*(1.d0-(max(emin0,etrans1)/etrans2)**(brad2+1.d0)) &
    /(brad2+1.d0)*etrans2
   endif
  else
   aw2=0.d0
  endif
  
  if(emin0.lt.etrans3)then
   gnorm3=etrans3*w_EBL_density(etrans3,zz)*sigics(e0,ame**2+etrans3*de)
   if(etrans3*de.gt.10.d0*ame**2)then
    gnorm3=max(gnorm3,w_EBL_density(etrans2,zz)*sigics(e0,ame**2 &
    +etrans2*de)*etrans2*exp(-(etrans3-etrans2)/akt/(1.d0+zz)))
    aw3=gnorm3*(exp((etrans3-max(emin0,etrans2))/akt/(1.d0+zz))-1.d0) &
    *akt*(1.d0+zz)
   else
    gnorm3=max(gnorm3,w_EBL_density(max(2.d0*emin0,etrans2),zz) &
    *sigics(e0,ame**2+max(2.d0*emin0,etrans2)*de)*etrans3 &
    *exp(-(etrans3-max(2.d0*emin0,etrans2))/akt/(1.d0+zz)))
    aw3=gnorm3*akt*(1.d0+zz)*(max(emin0,etrans2)/etrans3 &
    *exp((etrans3-max(emin0,etrans2))/akt/(1.d0+zz))-1.d0 &
    +akt*(1.d0+zz)/etrans3 &
    *(exp((etrans3-max(emin0,etrans2))/akt/(1.d0+zz))-1.d0))
   endif
  else
   aw3=0.d0
  endif
  
  gnorm4=eirmax*w_EBL_density(eirmax/2.d0,zz)*sigics(e0,ame**2+eirmax*de)
  if(etrans3*de.gt.10.d0*ame**2)then
   gnorm4=max(gnorm4,w_EBL_density(etrans3,zz) &
   *sigics(e0,ame**2+etrans3*de)*etrans3*(etrans3/eirmax)**3)
   aw4=gnorm4*((eirmax/max(emin0,etrans3))**2-1.d0)/2.d0*eirmax
  else
   gnorm4=max(gnorm4,w_EBL_density(etrans3,zz) &
   *sigics(e0,ame**2+etrans3*de)*etrans3*(etrans3/eirmax)**2)
   aw4=gnorm4*(eirmax/max(emin0,etrans3)-1.d0)*eirmax
  endif

1 rrr=psran()*(aw1+aw2+aw3+aw4)

! sample emin (= (sgam-m_e^2)/de) according to the 'proposal function';
! define 'rejection function' ( gb = f(emin) / f_proposal(emin) )
  if(rrr.lt.aw1)then
   emin=etrans1*(1.d0-psran()*(1.d0-(emin0/etrans1)**(brad1+1.d0))) &
   **(1.d0/(brad1+1.d0))
   gb0=gnorm1*(emin/emin0)**brad1
  elseif(rrr.lt.aw1+aw2)then
   if(etrans2*de.gt.10.d0*ame**2)then
    emin=etrans2*(1.d0-psran() &
    *(1.d0-(max(emin0,etrans1)/etrans2)**brad2))**(1.d0/brad2)
    gb0=gnorm2*(emin/etrans2)**(brad2-1.d0)
   else
    emin=etrans2*(1.d0-psran()*(1.d0-(max(emin0,etrans1)/etrans2) &
    **(1.d0+brad2)))**(1.d0/(1.d0+brad2))
    gb0=gnorm2*(emin/etrans2)**brad2
   endif

  elseif(rrr.lt.aw1+aw2+aw3)then
2  emin=etrans3-akt*(1.d0+zz)*dlog(1.d0-psran() &
   *(1.d0-exp((etrans3-max(emin0,etrans2))/akt/(1.d0+zz))))
   gb0=gnorm3*exp((etrans3-emin)/akt/(1.d0+zz))
   
   if(etrans3*de.lt.10.d0*ame**2)then
    if(psran().gt.emin/etrans3)goto 2
    gb0=gb0*emin/etrans3
   endif

  else
   if(etrans3*de.gt.10.d0*ame**2)then
    emin=eirmax/dsqrt(1.d0-psran()*(1.d0-(eirmax/max(emin0,etrans3))**2))
    gb0=gnorm4*(eirmax/emin)**3
   else
    emin=eirmax/(1.d0-psran()*(1.d0-eirmax/max(emin0,etrans3)))
    gb0=gnorm4*(eirmax/emin)**2
   endif
  endif
  
  gb0=gb0*gbnorm
  
  if(emin0/etrans3.gt..8d0) gb0=gb0  *5.d0
  
  if(model.eq.3.and.emin0/etrans3.gt..5d0) gb0=gb0*50.d0
    
  if(.not.(gb0.gt.0.d0.and.gb0.lt.1.d60))then
     write(*,*)'electron: gb0=0',gb0,gbnorm
     stop
     ierr=1
     return
  endif
  
  sgam=ame**2+emin*de                 ! c.m. energy for e(+-)-gamma interaction
  gb=w_EBL_density(emin,zz)*sigics(e0,sgam)*emin/gb0

  if (gb.gt.1.d0.and.nrenorm.eq.0) write(*,*)'sample_cmb(electron): gb='  &
  ,gb,nrej,nrenorm,emin0/etrans3
  
  if (psran().gt.gb) then       ! rejection
     nrej=nrej+1                ! total number of rejections for current sampling
     gbmax=max(gbmax,gb)        ! maximal value for rejection function
     if(nrej.gt.nrejmax)then    ! too many rejections
      if(gbmax.le.0.d0)then     ! wrong kinematics
       write(*,*)'electron: gbmax=0!!!'
       ierr=1
       return
      else
!       write(*,*)'nrej(e)>nrejmax'  &
!       ,nrej,emin0/etrans3,nrenorm,gbmax,e0/1.d9,emin0/etrans2  !,gb,gbnorm,gb0
       gbnorm=gbnorm*gbmax*2.d0 ! change normalization for the rejection function
       gbmax=0.d0
       nrenorm=nrenorm+1
       if(nrenorm.gt.100)stop
       nrej=0
      endif
     endif
     goto 1                     ! new try
  end if

end subroutine sample_electron
!=============================================================================!
!=============================================================================!
!     weighted background photon density Eq. (9)    (interpolation)           ! 
!=============================================================================!
double precision function w_EBL_density(emin,zz)
!-----------------------------------------------------------------------------
! input:
!        emin - minimal required energy for a background photon
!        zz - current redshift
!-----------------------------------------------------------------------------
!f2py intent(in) emin
!f2py intent(in) zz
!f2py intent(out) w_EBL_density
  use EBL_fit
  use xsec
  implicit none
  integer k,k1,jz,l1
  double precision emin,zz,dz,wk(3),wz(3),yl

  w_EBL_density=0.d0
  if (emin.ge.eirmax) return
      
  if (emin.le.erad1) then
     yl=log(emin*1.d15)/dlog(erad1*1.d15)*10.d0+1.d0
  else
     yl=log(emin/erad1)/dlog(eirmax/erad1)*40.d0+11.d0
  endif
  k=min(int(yl),48)
  k=max(k,1)
  if (k.eq.10) k=9
  wk(2)=yl-k
  wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
  wk(1)=1.d0-wk(2)+wk(3)
  wk(2)=wk(2)-2.d0*wk(3)
      
  dz=2.d0*zz+1.d0
  jz=min(9,int(dz))
  wz(2)=dz-jz
  wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
  wz(1)=1.d0-wz(2)+wz(3)
  wz(2)=wz(2)-2.d0*wz(3)

  do k1=1,3
     do l1=1,3
        w_EBL_density=w_EBL_density+denc(k+k1-1,jz+l1-1)*wk(k1)*wz(l1)
     end do
  end do
  w_EBL_density=exp(w_EBL_density-(1.d0-brad1)*dlog(emin))

end function w_EBL_density
!=============================================================================!
!=============================================================================!
!  interaction length for interactions on EBL photons    (interpolation)      !
!=============================================================================!
double precision function int_length(e0,x0,zz,icq)
!-----------------------------------------------------------------------------
! calls: rate_bb, psran
! input:
!        e0 - particle energy;
!        x0 - particle distance from the source;
!        zz - particle redshift position;
!        icq - particle type (0 - photon, +/-1 - electron/positron)
!-----------------------------------------------------------------------------
!f2py intent(in) e0
!f2py intent(in) x0
!f2py intent(in) zz
!f2py intent(in) icq
!f2py intent(out) int_length
  use internal, only : rcmb
  implicit none
  double precision e0,x0,zz,sig
  double precision rate_bb,psran
  integer icq
      
  sig = rate_bb(e0,zz,icq)           ! inverse mean free pass
  if (sig.le.0.d0) then              ! particle energy below threshold
     int_length=1.1d0*rcmb           ! -> no interaction
  else
     int_length=x0-log(psran())/sig  ! interaction lengt 
  endif
  
end function int_length
!=============================================================================!
!=============================================================================!
!        pair production x-section                                            !
!=============================================================================!
double precision function sigpair(sgam)
!-----------------------------------------------------------------------------
! input:
!        sgam - c.m. energy for gamma-gamma interaction
!-----------------------------------------------------------------------------
!f2py intent(in) sgam
!f2py intent(out) sigpair
  use constants, only : ame,sigtmp
  implicit none
  double precision sgam,bet
      
  bet=1.d0-4.d0*ame**2/sgam
  if (bet.le.0.d0) then
     sigpair=0.d0
  else
     bet=sqrt(max(0.d0,bet))
     sigpair=sigtmp*3.d0/4.d0*ame**2/sgam*((3.d0-bet**4) &
 & *log((1.d0+bet)/(1.d0-bet))-2.d0*bet*(2.d0-bet**2))
  end if

end function sigpair
!=============================================================================!
!=============================================================================!
!        energy partition for pair production                                 !
!=============================================================================!
double precision function zpair(sgam)
!-----------------------------------------------------------------------------
! calls: psran
! input:
!        sgam - c.m. energy for gamma-gamma interaction
!-----------------------------------------------------------------------------
!f2py intent(in) sgam
!f2py intent(out) zpair
  use constants, only : ame
  use internal, only : debug
  implicit none
  double precision sgam,bet,zmin,z,gb
  double precision psran

  bet=sqrt(max(0.d0,1.d0-4.d0*ame**2/sgam))
  zmin=(1.d0-bet)/2.d0
  do 
     z=.5d0*(2.d0*zmin)**psran()
     gb=(z**2/(1.d0-z)+1.d0-z+(1.d0-bet**2)/(1.d0-z)-(1.d0-bet**2)**2 &
       &/4.d0/z/(1.d0-z)**2)/(1.d0+2.d0*bet**2*(1.d0-bet**2))
     if (debug>0.and.gb.gt.1.d0) write(*,*)'zpair: gb=',gb
     if (psran()<gb) exit
  end do
  if (psran()>0.5d0) z=1.d0-z
  zpair=z     

end function zpair
!=============================================================================!
!=============================================================================!
!        (inverse) Compton x-section                                          !
!=============================================================================!
double precision function sigics(e0,sgam)
!-----------------------------------------------------------------------------
! input:
!        e0 - electron/positron energy;
!        sgam - c.m. energy for e(+-)-gamma interaction
!-----------------------------------------------------------------------------
!f2py intent(in) sgam
!f2py intent(in) e0
!f2py intent(out) sigics
  use constants, only : ame,sigtmp
  use user_variables, only : ethr
  implicit none
  double precision e0,sgam,zmax,zm,t
      
  zmax=1.d0-ethr/e0
  zm=ame**2/sgam
  zmax=max(zmax,zm)
  if (1.d0-zm.le.3.d-3) then              ! use series expansion near threshold
     t=min(1.d0,(zmax-zm)/(1.d0-zm))
     sigics=.75d0*sigtmp*t*(4.d0*t*(t*(.5d0/zmax-(1.d0+zm)/3.d0 &
  &  /zm**2)-min(1.d0,(1.d0-zmax)/(1.d0-zm))*zm/zmax/2.d0) &
  &  +(zmax+zm)*zm/2.d0+1.d0-(zmax-zm)/zm/2.d0 &
  &  +(zmax-zm)**2/zm**2/3.d0)
  else
     sigics=.75d0*sigtmp*zm*min(1.d0,(zmax-zm)/(1.d0-zm)) &
  &  *(dlog(zmax/zm)/(zmax-zm)*(1.d0-4.d0*zm*(1.d0+zm)/(1.d0-zm)**2) &
  &  +4.d0*(zm/zmax+zm)/(1.d0-zm)**2+(zmax+zm)/2.d0)
  end if
  if(sigics.le.0.d0)then
  ! write(*,*)'sigics:e0,sgam,zm,zmax',e0,sgam,zm,zmax,sigics
   sigics=0
  endif
end function sigics
!=============================================================================!
!=============================================================================!
!        energy partition for inverse Compton    (E-fraction taken by e+-)    !
!=============================================================================!
double precision function zics(e0,sgam)
!-----------------------------------------------------------------------------
! calls: psran
! input:
!        sgam - c.m. energy for e(+-)-gamma interaction
!-----------------------------------------------------------------------------
!f2py intent(in) sgam
!f2py intent(out) zics
  use constants, only : ame
  use user_variables, only : ethr
  use internal, only : debug
  implicit none
  double precision e0,sgam,zmin,zmax,z,gb
  double precision psran

  zmax=1.d0-ethr/e0
  zmin=ame**2/sgam

  if (zmin.ge.zmax) then
     if (debug>0) call error('zmin>zmax in zics',1)
     zics=zmin
  else
     do
        z=zmin*(zmax/zmin)**psran()
        gb=(1.d0+z*z)/2.d0-2.d0*zmin/z*(z-zmin)*(1.d0-z)/(1.d0-zmin)**2
!!!        if (debug>0.and.gb.gt.1.d0) write(*,*)'zics: gb=',gb,z,zmin,zmax
        if (psran()<gb) exit
     end do
     zics=z
  endif
  
end function zics
!=============================================================================!
!=============================================================================!
!        cross section * energy loss (under threshold) for ICS
!=============================================================================!
double precision function zsigics(e0,sgam) 
!-----------------------------------------------------------------------------
! calls: error
! input:
!        e0 - electron/positron energy;
!        sgam - c.m. energy for e(+-)-gamma interaction
!-----------------------------------------------------------------------------
!f2py intent(in) sgam
!f2py intent(in) e0
!f2py intent(out) zsigics
  use constants, only : ame,sigtmp
  use user_variables, only : ethr
  implicit none
  double precision e0,sgam,zmax,zm
      
  zm=ame**2/sgam
  zmax=max(1.d0-ethr/e0,zm)
      
  if(zm.lt..5d0)then
     zsigics=.75d0*sigtmp*zm*(1.d0-zmax)/(1.d0-zm)                   &
  &  *((-dlog(zmax)/(1.d0-zmax)-1.d0)*(1.d0-4.d0*zm*(1.d0+2.d0*zm)   &
  &  /(1.d0-zm)**2)+(1.d0-zmax)*(1.d0+2.d0*zmax)/6.d0                &
  &  +2.d0*zm*(1.d0-zmax)/(1.d0-zm)**2*(1.d0+2.d0*zm/zmax))
  else
     zsigics=.75d0*sigtmp*zm*(1.d0-zmax)**2/(1.d0-zm)                &
  &  *(1.d0+.25d0*(1.d0-zmax)**2+4.d0*zm*(1.d0-zmax)/(1.d0-zm)**2    &
  &  *(zm/zmax-(1.d0+2.d0*zm)/3.d0*(1.d0+.75d0*(1.d0-zmax))))
  endif
     
  if(zsigics.lt.0.d0)then
     call error('zsigics<0',1)
     zsigics=0.d0
  endif

end function zsigics
!=============================================================================!
!=============================================================================!
!        relative energy loss (per cm) (due to under-threshold ICS photons)   !
!=============================================================================!
double precision function zloss(e0,zz)
!-----------------------------------------------------------------------------
! input:
!        e0 - electron/positron energy;
!        zz - current redshift
!-----------------------------------------------------------------------------
!f2py intent(in) zz 
!f2py intent(in) e0
!f2py intent(out) zloss
  use constants
  use EBL_fit
  use user_variables, only : ethr,egmax
  use xsec
  implicit none
  integer jz,k,k1,l1
  double precision e0,zz,yl,dz,wk(3),wz(3)

  zloss=0.d0
  yl=dlog(e0/ame)/dlog(egmax/ame)*50.d0+1.d0
  k=min(int(yl),49)
  k=max(k,2)
  wk(2)=yl-k
  wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
  wk(1)=1.d0-wk(2)+wk(3)
  wk(2)=wk(2)-2.d0*wk(3)

  dz=2.d0*zz+1.d0
  jz=min(9,int(dz))
  wz(2)=dz-jz
  wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
  wz(1)=1.d0-wz(2)+wz(3)
  wz(2)=wz(2)-2.d0*wz(3)

  do k1=1,3
     do l1=1,3
        zloss=zloss+zsigc(k+k1-1,jz+l1-1)*wk(k1)*wz(l1)
     enddo
  enddo
  zloss=exp(zloss)

end function zloss
!=============================================================================!
!=============================================================================!
!       interpolates interaction rate/cm on photon background                 !
!=============================================================================!
double precision function rate_bb(e0,zz,icq)
!-----------------------------------------------------------------------------
! input:
!        e0 - particle energy;
!        zz - particle redshift position;
!        icq - particle type (0 - photon, +/-1 - electron/positron)
!-----------------------------------------------------------------------------
!f2py intent(in) zz 
!f2py intent(in) e0
!f2py intent(in) icq
!f2py intent(out) rate_bb
  use constants
  use EBL_fit
  use user_variables, only : ethr,egmax
  use xsec
  implicit none
  integer icq,ica,jz,k,k1,l1
  double precision e0,zz,emin,yl,dz,wk(3),wz(3)

  rate_bb=0.d0
  if (icq.eq.0) then
     emin=ame**2/eirmax
  else
     emin=ethr+.5d0*ame**2/eirmax/(1.d0+dsqrt(1.d0+ame**2/ethr**2/eirmax*(ethr-eirmax)))
  end if
  if (e0.le.emin) return
      
  ica=iabs(icq)+1
  yl=log(e0/emin)/log(egmax/emin)*50.d0+1.d0
  k=min(int(yl),49)
  k=max(k,2)
  wk(2)=yl-k
  wk(3)=wk(2)*(wk(2)-1.d0)*.5d0
  wk(1)=1.d0-wk(2)+wk(3)
  wk(2)=wk(2)-2.d0*wk(3)

  dz=2.d0*zz+1.d0
  jz=min(9,int(dz))
  wz(2)=dz-jz
  wz(3)=wz(2)*(wz(2)-1.d0)*.5d0
  wz(1)=1.d0-wz(2)+wz(3)
  wz(2)=wz(2)-2.d0*wz(3)

  do k1=1,3
     do l1=1,3
        rate_bb=rate_bb+sigc(k+k1-1,jz+l1-1,ica)*wk(k1)*wz(l1)
     enddo
  enddo
  rate_bb=exp(rate_bb)
  
end function rate_bb
!=============================================================================!
!=============================================================================!
!        electron synchrotron energy loss (eV/cm)                     !
!=============================================================================!
double precision function eloss(e0,begmf)
!-----------------------------------------------------------------------------
! input:
!        e0 - electron/positron energy;
!        begmf - strength of (transverse) extragalactic B-field
!-----------------------------------------------------------------------------
!f2py intent(in) begmf
!f2py intent(in) e0
!f2py intent(out) eloss
  use constants
  implicit none
  double precision e0,begmf,chi

  chi=sqrt(max(0.d0,(e0/ame)**2-1.d0))*begmf/bcr
  eloss = chi**2/(1.d0+4.8d0*(1.d0+chi)*log(1.d0+1.7d0*chi) +  &
           3.44d0*chi**2)**(2.d0/3.d0)*ame**2/137.d0/1.5d0*1.d7/197.d0 

end function eloss
!=============================================================================!
!=============================================================================!
!        electron/positron deflection angle (rad/cm) by EGMF                  !
!=============================================================================!
double precision function themf(e0,begmf)
!-----------------------------------------------------------------------------
! input:
!        e0 - electron/positron energy;
!        begmf - strength of (transverse) extragalactic B-field
!-----------------------------------------------------------------------------
!f2py intent(in) begmf
!f2py intent(in) e0
!f2py intent(out) themf
  implicit none
  double precision e0,begmf
  themf=294.d0*begmf/e0
end function themf
!=============================================================================!
!=============================================================================!
