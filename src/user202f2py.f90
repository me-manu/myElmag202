!=============================================================================!
!=============================================================================!
module user_result
  implicit none
  save
  double precision ::        &
     log10t_min=-3.00d0,      &
     log10t_max=9.00d0,      &
     log10theta_min=-8.00d0, &
     log10theta_max=1.30d0
  integer ::    &
     n_bine=81,  &
     n_binth=93, &
     n_bint=48
  double precision, allocatable, dimension (:,:) :: hist_ge  ! diffuse spectrum
  double precision, allocatable, dimension (:,:) :: hist_eth ! spectrum
  double precision, allocatable, dimension (:,:) :: hist_et  ! spectrum
  double precision, allocatable, dimension (:,:,:) :: hist_etht ! spectrum
  integer n_reg
end module user_result
!=============================================================================!
!=============================================================================!
! Variables for the cascade determined by the user                            !
! In the standard ELMAG version these are all parameters, which leads         !
! to the fact that the parameters cannot be changed from python               !
!=============================================================================!
module user_variables
  implicit none
  save
  integer :: model=6
                                   ! 1: Kneiske Dole best fit, 2: lower limit,
                                   ! 3: Franceschini arXiv:0805.1841
                                   ! 4: Model C Finke et al. arXiv:0905.1115
                                   ! 5: Gilmore et al. arXiv:1104.0671
                                   ! 6: Dominguez et al. arXiv:

  integer :: nmax=1*10**3 ! number of simulated particles --Needed??

  double precision ::        &  
  ethr=1.00d8, &
  egmax=15.d12, &
  ir_rescale=1.00d0, &
  cohlnth=3.086d21*1.d3*1.00d0, &
  th_jet=6.00d0, &
  a_smp=0.00d0
  
  double precision :: igmf=1.00d-16

  ! gamma spectrum
  double precision ::        &  
  emin=1.00d9, &
  ebreak=1.50d13, &
  gam1=-1.70d0, &
  gam2=-1.70d0
  character*400 :: tablefile_n = 'Tables/n_Dom.dat' ! file for photon density file
  character*400 :: tablefile_z = 'Tables/z-IR_Dom.dat' ! file for redshifts if required
end module user_variables

!=============================================================================!
!=============================================================================!
!=============================================================================!
! sample initial photon/electron from the input spectra                       !
!=============================================================================!
subroutine initial_particle(e0,weight)
!-----------------------------------------------------------------------------
! input:
!        e0 - particle energy;
!        weight - initial particle weight
!-----------------------------------------------------------------------------
!f2py intent(out) e0 
!f2py intent(out) weight
  use user_variables
  implicit none
  double precision e0,weight
  double precision psran
      
  e0=emin*(egmax/emin)**psran()                      !energy: uniform in ln E 
  if (e0<ebreak) then
     weight=(e0/ebreak)**(gam1+1.d0)*log(egmax/emin)
     !weight=(e0/ebreak)**(gam1)*log(egmax/emin)
  else
     weight=(e0/ebreak)**(gam2+1.d0)*log(egmax/emin) 
     !weight=(e0/ebreak)**(gam2)*log(egmax/emin) 
  endif

end subroutine initial_particle
!=============================================================================!
!=============================================================================!
!           strength of EGMF/Gauss                                            !
!=============================================================================!
double precision function bemf(r)
!-----------------------------------------------------------------------------
! input:
!        r - distance to Earth/cm
!-----------------------------------------------------------------------------
!f2py intent(in) r
!f2py intent(out) bemf
  use user_variables, only : igmf
  implicit none
  double precision r
  bemf = igmf
end function bemf
!=============================================================================!
!=============================================================================!
subroutine registerresult(e0,theta,dt,weight,icq)
!f2py intent(in) e0
!f2py intent(in) theta 
!f2py intent(in) dt
!f2py intent(in) weight
!f2py intent(in) icq 
  use user_variables
  use user_result
  implicit none
  integer icq,i,j,k,tb
  double precision e0,weight,theta,dt
  double precision add
  REAL, DIMENSION(n_binth + 1) :: theta_bounds
  REAL, DIMENSION(n_bint + 1) :: t_bounds
  theta_bounds = (/ (log10theta_min + I * (log10theta_max - (log10theta_min)) / n_binth,  I = 0, n_binth) /)
  t_bounds = (/ (log10t_min + I * (log10t_max - (log10t_min)) / n_bint,  I = 0, n_bint) /)

  !add = weight*e0/log(egmax/ethr)*(n_bine-1)
  add = weight/log(egmax/ethr)*(n_bine-1)
  ! Weight definition: weight=(e0/ebreak)**(gam1+1.d0)*log(egmax/emin)
! diffuse energy spectrum:
  i=min(n_bine,int(log(e0/ethr)/log(egmax/ethr)*(n_bine-1))+1)
  i=max(i,1)
  hist_ge(i,abs(icq)+1)=hist_ge(i,abs(icq)+1)+add


  if (icq.ne.0) return                                 ! forget electrons


  if (log10(theta)<theta_bounds(1)) then
      j = 1
  elseif (log10(theta)>theta_bounds(n_binth)) then
      j = n_binth
  else
      do j = 2, size(theta_bounds)
      if (log10(theta) >= theta_bounds(j - 1) .AND. log10(theta) < theta_bounds(j)) then
      exit
      endif
      end do
  endif

  !!write(*,*) 'log10 theta',log10(theta),j
  hist_eth(i,j)=hist_eth(i,j)+add
  !write(*,*) 'register e0,theta,dt,weight,icq,i,j', e0,log10(theta),log10(dt),weight,icq,i,j
  !write(*,*) 'e0,theta,dt', e0,log10(theta),log10(dt)
  !write(*,*) e0,log10(theta),log10(dt)

  if (log10(dt)< t_bounds(1)) then
      k = 1
  elseif (log10(dt)>t_bounds(n_bint)) then
      k = n_bint
  else
      do k = 2, size(t_bounds)
      if (log10(dt) >= t_bounds(k - 1) .AND. log10(dt) < t_bounds(k)) then
      exit
      endif
      end do
  endif
  !!write(*,*) t_bounds
  !!write(*,*) 'dt',log10(dt), k
  hist_et(i,k)=hist_et(i,k)+add
  hist_etht(i,j,k)=hist_etht(i,j,k)+add

! angular profile:
  n_reg = n_reg+1

end subroutine registerresult
!=============================================================================!
