!=============================================================================!
!=============================================================================!
module stack
  use, intrinsic :: iso_fortran_env
  implicit none
  save
  integer, parameter ::  n_max = 100             ! max. number of secondaries
  integer :: jcmb

  type one_event
     !integer(int8) :: icq
     integer :: icq
     double precision :: en,x,z,t,w,the1,the2,xxc,xx,dt
     !real(real64) :: en,x,z,t,w,the1,the2,xxc,xx,dt
  end type one_event 
  integer, parameter :: event_size = 8 + 10 * 64
  !f2py integer*1, dimension(event_size) :: act
  !f2py integer*1, dimension(event_size , n_max) :: event
  type (one_event) act
  type (one_event) event(n_max)

end module stack
!=============================================================================!
!=============================================================================!
module EBL_fit  
  implicit none
  save
  integer :: n_EBL1, n_EBL2

  double precision, parameter ::  akt=2.348d-4,eirmin=2.5d-3,eirmax=12.d0,  &
   erad1=4.136d-9,erad2=1.05d-6,arad1=1.039d15,arad2=4.408d-7,brad1=0.6d0, &
   brad2=-1.65d0,crad=1.047d18,drad=1.85d0

  double precision, allocatable, dimension(:) :: z_ir,eir
  double precision, allocatable, dimension(:,:) :: anir
 
end module EBL_fit
!=============================================================================!
!=============================================================================!
module xsec        
  implicit none
  save
  double precision sigc(51,11,2),denc(51,11),zsigc(51,11)
end module xsec
!=============================================================================!
!=============================================================================!
module constants
  implicit none
  save
  double precision, parameter ::             & 
    pi=3.1415926536d0,                       &
    two_pi=2.d0*pi,                          &
    degree_rad = pi/180.d0,                  &     
    rad_degree = 180.d0/pi

  double precision, parameter ::             & 
    bcr=4.14d13,                             & !critical magnetic field, Gauss
    ame=5.11d5,                              & !electron mass
    sigtmp=6.6525d-25,                       & !Thomson x-section, cm^2
    rthomp=9.17d11,                          & !Thomson length, cm
    cy=9.46d17                                 !speed of light, cm/year

end module constants
!=============================================================================!
!=============================================================================!
module cosmology
  implicit none
  save
  double precision, parameter ::          & 
  T_CMB = 2.35d-4,                        & ! CMB temperature in eV
  Omega_m = 0.3d0,                        & ! matter
  Omega_v = 0.7d0,                        & ! vacuum energy
  h = 0.7d0,                              & ! Hubble constant 
  H_0 = 1.d0/3.0856d17 * h,               &
  Mpc_s = 1.d14,                          & ! Mpc in s
  R_z3=6753.5d0                             ! comoving radial distance at z=3

  double precision :: t_0

end module cosmology
!=============================================================================!
!=============================================================================!
module internal
  implicit none
  save
  integer, parameter :: K4B=selected_int_kind(9)
  integer(K4B) :: iseed
  double precision bal,rcmb
  integer, parameter ::                     & 
    debug = 1                             ! 1,2 performs tests, prints infos
end module internal
!=============================================================================!
!=============================================================================!
