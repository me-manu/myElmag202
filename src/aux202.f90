!.............................................................................
!                         error handling                                     .
!.............................................................................
subroutine error(estring,s)
!f2py intent(in) estring
!f2py intent(in) s
  character (len=*), intent(in) :: estring
  integer, intent(in) :: s
  integer, save :: n_warn  
  integer :: n_warn_max = 100

  if (s==1.or.s==11) then                   ! warning message
     write(*,*) 
     write(*,*) 'Warning:'
     write(*,*) estring
     write(99,*) estring
     if (s==1) then                         ! warning 
        n_warn = n_warn+1
        if (n_warn>n_warn_max) then
           write(*,*)
           write(*,*) 'more than',n_warn_max,' warnings!'
           write(*,*)
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*) '!  ELMAG 1.01 stops program excecution  !'
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           stop
        endif
     endif
  endif
  if (s==0) then                    ! error
     write(*,*)
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) '!   a serious error:                    !'
     write(99,*)'!   a serious error:                    !'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*)
     write(*,*) estring
     write(99,*)estring
     write(*,*)
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) '!  ELMAG 1.01 stops program excecution  !'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     stop
  endif

end subroutine error

!..............................................................................
!          random number generator from Numerical Recipes (Fortran90)         .
!..............................................................................
function ran0()
!f2py intent(out) ran0
  use internal, only : iseed
  implicit none
  integer, parameter :: K4B=selected_int_kind(9)
  double precision ran0
  integer(K4B),parameter :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  real, save :: am
  integer(K4B), save :: ix=-1, iy=-1,k

  if (iseed <= 0 .or. iy < 0) then
     am = nearest(1.e0,-1.e0)/IM
     iy = ior(ieor(888889999,abs(iseed)),1)
     ix = ieor(777755555,abs(iseed))
     iseed = abs(iseed)+1
  end if
  ix = ieor(ix,ishft(ix,13))
  ix = ieor(ix,ishft(ix,-17))
  ix = ieor(ix,ishft(ix,5))
  k = iy/IQ
  iy = IA*(iy-k*IQ)-IR*k
  if (iy<0) iy = iy+IM
  ran0 = am*ior(iand(IM,ieor(ix,iy)),1)
  !write(*,*) ran0, iseed
!..............................................................................
end function ran0
double precision function psran()
!f2py intent(out) psran
  implicit none
  double precision ran0
  psran=ran0()
end function psran
!..............................................................................
!..............................................................................
! Performs polynomial interpolation. The first term F(NN) contains the  values.
! of the function to interpolate, A(NN) contains the corresponding x-values,  .
! x is the point at which we want to evaluate the function and the last 
! parameter MM corresponds to the order of the polynomial 
! interpolation 1<MM<9 (set MM=3 if you don't know). This interpolating 
! routine requires real numbers in singol precision and an ordering
! AA(i)<AA(i+1) as input.
!..............................................................................
double precision function DIVDIF(F,A,NN,X,MM)
  implicit none
  integer n,m,nn,mm,mplus,i,ip,ix,iy,isub,j,l,mid,npts,mmax
  double precision x,A(NN),F(NN),T(20),D(20),sum
  LOGICAL EXTRA
!f2py intent(in) F
!f2py intent(in) A
!f2py intent(in) NN
!f2py intent(in) X
!f2py intent(in) MM
!f2py intent(out) DIVDIF
!f2py depend(NN) F
!  LOGICAL MFLAG,RFLAG
  DATA MMAX/10/
!
!  TABULAR INTERPOLATION USING SYMMETRICALLY PLACED ARGUMENT POINTS.
!
!  START.  FIND SUBSCRIPT IX OF X IN ARRAY A.
  IF( (NN.LT.2) .OR. (MM.LT.1) ) GO TO 20
  N=NN
  M=MM
  MPLUS=M+1
  IX=0
  IY=N+1
 !     (SEARCH INCREASING ARGUMENTS.)
1 MID=(IX+IY)/2
  IF(X.GE.A(MID)) GO TO 2
  IY=MID
  GO TO 3
!        (IF TRUE.)
2 IX=MID
3 IF(IY-IX.GT.1) GO TO 1
  GO TO 7
!  COPY REORDERED INTERPOLATION POINTS INTO (T(I),D(I)), SETTING
!  *EXTRA* TO TRUE IF M+2 POINTS TO BE USED.
7 NPTS=M+2-MOD(M,2)
  IP=0
  L=0
  GO TO 9
8 L=-L
  IF(L.GE.0) L=L+1
9 ISUB=IX+L
  IF((1.LE.ISUB).AND.(ISUB.LE.N)) GO TO 10
!        (SKIP POINT.)
  NPTS=MPLUS
  GO TO 11
!        (INSERT POINT.)
10 IP=IP+1
  T(IP)=A(ISUB)
  D(IP)=F(ISUB)
11 IF(IP.LT.NPTS) GO TO 8
  EXTRA=NPTS.NE.MPLUS
!
!  REPLACE D BY THE LEADING DIAGONAL OF A DIVIDED-DIFFERENCE TABLE, SUP-
!  PLEMENTED BY AN EXTRA LINE IF *EXTRA* IS TRUE.
  DO L=1,M
     IF(.NOT.EXTRA) GO TO 12
     ISUB=MPLUS-L
     D(M+2)=(D(M+2)-D(M))/(T(M+2)-T(ISUB))
12   I=MPLUS
     DO J=L,M
        ISUB=I-L
        D(I)=(D(I)-D(I-1))/(T(I)-T(ISUB))
        I=I-1
     end DO
  end DO        
!
!  EVALUATE THE NEWTON INTERPOLATION FORMULA AT X, AVERAGING TWO VALUES
!  OF LAST DIFFERENCE IF *EXTRA* IS TRUE.
  SUM=D(MPLUS)
  IF(EXTRA) SUM=0.5*(SUM+D(M+2))
  J=M
  DO L=1,M
     SUM=D(J)+(X-T(J))*SUM
     J=J-1
  end DO
  DIVDIF=SUM
  RETURN

20 IF(MM.LT.1) WRITE(*,101) MM
  IF(NN.LT.2) WRITE(*,102) NN

!            DIVDIF=999999999999999999.
  DIVDIF=999999999.e9

101 FORMAT( 7X, 'FUNCTION DIVDIF ... M =',I6,' IS LESS THAN 1')
102 FORMAT( 7X, 'FUNCTION DIVDIF ... N =',I6,' IS LESS THAN 2')

end function DIVDIF


     
     

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!            injection time t as function of source distance r                .
double precision function t_r(r)
!f2py intent(in) r
!f2py intent(out) t_r
  implicit none
  integer n,n1
  parameter(n=310,n1=311)
  double precision r,t_i(0:n),r_i(0:n)
  double precision divdif
  common /data_div3/ t_i,r_i
	
  t_r = DIVDIF(t_i,r_i,n1,r,3) 

end function t_r

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
!               source distance r as function of redshift z 
!------------------------------------------------------------------------------!
double precision function r_z(z)
!f2py intent(in) z
!f2py intent(out) r_z
  implicit none
  integer n,n1
  parameter(n=310,n1=311)
  double precision z,z_i(0:n),r_i(0:n)
  double precision divdif
  common /data_div4/ z_i,r_i
	
  r_z = DIVDIF(r_i,z_i,n,z,3) 

end function r_z

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!             redshift z  as function of source distance r 

double precision function z_r(r)
!f2py intent(in) r
!f2py intent(out) z_r
  implicit none
  integer n,n1
  parameter(n=310,n1=311)
  double precision r,r_i(0:n),z_i(0:n)
  double precision divdif
  common /data_div4/ z_i,r_i
        
  z_r = DIVDIF(z_i,r_i,n1,r,3) 

end function z_r


!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!            age of a flat Universe as function of redshift z                 .
!------------------------------------------------------------------------------!
double precision function t_z(z)
!f2py intent(in) z
!f2py intent(out) t_z
  implicit none
  integer n
  parameter(n=310)
  double precision t_i(0:n),z_i(0:n),xp,yp
  double precision divdif
  double precision z
  common /data_div1/ t_i,z_i
  
  xp = z
  yp = DIVDIF(t_i,z_i,n,xp,3) 
  t_z = yp           

end function t_z

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------------!
!            redshift z of a flat Universe as function of age t/t_0           .
!------------------------------------------------------------------------------!
double precision function z_t(t)
!f2py intent(in) t
!f2py intent(out) z_t
  implicit none
  integer n,n1
  parameter(n=310,n1=311)
  double precision t,t_i(0:n),z_i(0:n)
  double precision divdif
  common /data_div2/ t_i,z_i

  z_t = DIVDIF(z_i,t_i,n1,t,3) 

end function z_t

