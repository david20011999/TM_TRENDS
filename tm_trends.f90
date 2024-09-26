module vers
  character(len=50):: version='TM - 02 August 2011'
  contains
    subroutine print_version(start)
    implicit none
    logical:: start
    
    
    print *,'-----------------'
    print *,trim(version)
    print *,'by A Legarra, L Varona, E Lopez de Maturana'
    print *,'-----------------'
    if(start) print *,'started: '
    call printtime2
    end subroutine

      subroutine printtime2
      INTEGER  DATE_TIME (8)
      CHARACTER(LEN = 12) REAL_CLOCK (3)

      CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &
                      REAL_CLOCK (3), DATE_TIME)

      print *,'date: ',real_clock(1)(7:8),'/',real_clock(1)(5:6),&
      '/',real_clock(1)(1:4)
      print *,'time: ',&
      real_clock(2)(1:2),':',real_clock(2)(3:4),':',real_clock(2)(5:6)

      end subroutine

  
end module vers

module kinds
! From Ignacy Misztal BLUPF90 distribution 
  integer, parameter :: single = SELECTED_REAL_KIND( 6, 37 )
  integer, parameter :: double = SELECTED_REAL_KIND( 15, 307 )
  integer, parameter :: extended = SELECTED_REAL_KIND( 18, 4931 )
  integer, parameter :: r4 = SELECTED_REAL_KIND( 6, 37 )
  integer, parameter :: r8 = SELECTED_REAL_KIND( 15, 307 )
  integer, parameter :: r16 = SELECTED_REAL_KIND( 18, 4931 )

  ! current precison for hash storage
  integer, parameter :: rh=r8
end module kinds

module modelo
  implicit none
      integer :: &
      maxrec=5000000,& !this is the number of non-zero elements in the MME
     iguard=1000000000 !never write continuation file
! total number of iterations : imue*icad
! burn-in: lap*icad
! thin interval or lag: icad
! number of samples used in estimating means and sd: (imue-lap)
! This is later read from par file


end module modelo

module listaligada
  use kinds
  USE modelo
  implicit none
  real(r8),allocatable:: zhz(:)
  integer,allocatable:: ifirst(:),ivcol(:),inext(:)
  integer :: nplace
 
  contains

!=======================================================================
      SUBROUTINE LINKS(IROW,ICOL,D)

!     CONSTRUCCION DE LISTAS LIGADAS POR FILAS
!     CON ORDEN CRECIENTE POR COLUMNAS
!     NO SE REQUIERE NINGUN ORDEN EN LOS DATOS


	  double precision:: D
	  integer:: ipre,iplace,irow,icol

      IPRE=0
      IPLACE=IFIRST(IROW)
4     IF(IPLACE.GT.0)THEN
        IF(IVCOL(IPLACE).GE.ICOL)THEN
          IF(IVCOL(IPLACE).EQ.ICOL)THEN
            ZHZ(IPLACE)=ZHZ(IPLACE)+D
            RETURN
          ELSE
            NPLACE=NPLACE+1
            IF(NPLACE.GT.MAXREC)THEN
              PRINT *,'There is no space in the linked list'
	      print *,'increasing linked list by 50%'
	      call increase_linked_list()
            END IF
            IF(IPRE.EQ.0)THEN
              INEXT(NPLACE)=IFIRST(IROW)
              IFIRST(IROW)=NPLACE
            ELSE
              INEXT(NPLACE)=INEXT(IPRE)
              INEXT(IPRE)=NPLACE
            END IF
            ZHZ(NPLACE)=D
            IVCOL(NPLACE)=ICOL
            RETURN
          END IF
        ELSE
          IPRE=IPLACE
          IPLACE=INEXT(IPLACE)
          GOTO 4
        END IF
      ELSE
        NPLACE=NPLACE+1
        IF(NPLACE.GT.MAXREC)THEN
              PRINT *,'There is no space in the linked list'
	      print *,'increasing linked list by 50%'
	      call increase_linked_list()
        END IF
        IF(IFIRST(IROW).GT.0)THEN
          INEXT(IPRE)=NPLACE
        ELSE
          IFIRST(IROW)=NPLACE
        END IF
        ZHZ(NPLACE)=D
        IVCOL(NPLACE)=ICOL
        INEXT(NPLACE)=0
      END IF
      END subroutine

!      -------------------------------------------------------------

 
  subroutine increase_linked_list()
  ! increase linked list structure by 50%
  ! simple because linked list is unordered
  ! AL 6/7/07
  implicit none
  integer,allocatable:: ivcol_temp(:),inext_temp(:)
  real(r8), allocatable:: zhz_temp(:)
  integer:: newmaxrec
    
  ! create temp storage
  allocate(ivcol_temp(size(ivcol)),&
           inext_temp(size(inext)),&
           zhz_temp(size(zhz))   )

  ivcol_temp=ivcol
  inext_temp=inext
  zhz_temp=zhz
         
  newmaxrec=floor(1.5*maxrec) !new size
  print *,'changing maxrec from: ', maxrec, 'to: ',newmaxrec
  deallocate(ivcol,&
             inext,&
        zhz)
  allocate(ivcol(0:newmaxrec), &
           inext(0:newmaxrec),&
           zhz(newmaxrec))	 
  ivcol=0
  inext=0
  zhz=0d0
  ivcol(0:maxrec)=ivcol_temp
  inext(0:maxrec)=inext_temp
  zhz(1:maxrec)=zhz_temp
  deallocate(ivcol_temp,&
             inext_temp,&
             zhz_temp)
  maxrec=newmaxrec	     
  
  end subroutine

end module listaligada

MODULE Ecuyer_random
! L'Ecuyer's 1996 random number generator.
! Fortran version by Alan.Miller @ vic.cmis.csiro.au
! N.B. This version is compatible with Lahey's ELF90
! http://www.ozemail.com.au/~milleraj
! Latest revision - 30 March 1999

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307) !=real(8) con DVF

! These are unsigned integers in the C version
INTEGER, SAVE :: s1 = 1234, s2 = -4567, s3 = 7890

CONTAINS

SUBROUTINE init_seeds(i1, i2, i3)
IMPLICIT NONE

INTEGER, INTENT(IN) :: i1, i2, i3

s1 = i1
s2 = i2
s3 = i3
IF (IAND(s1,-2) == 0) s1 = i1 - 1023
IF (IAND(s2,-8) == 0) s2 = i2 - 1023
IF (IAND(s3,-16) == 0) s3 = i3 - 1023

RETURN
END SUBROUTINE init_seeds


FUNCTION taus88() RESULT(random_numb)
! Generates a random number between 0 and 1.  Translated from C function in:
! Reference:
! L'Ecuyer, P. (1996) `Maximally equidistributed combined Tausworthe
! generators', Math. of Comput., 65, 203-213.

! The cycle length is claimed to be about 2^(88) or about 3E+26.
! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).

IMPLICIT NONE
REAL (dp) :: random_numb

INTEGER   :: b

! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
!      to the left if j > 0, otherwise to the right.

b  = ISHFT( IEOR( ISHFT(s1,13), s1), -19)
s1 = IEOR( ISHFT( IAND(s1,-2), 12), b)
b  = ISHFT( IEOR( ISHFT(s2,2), s2), -25)
s2 = IEOR( ISHFT( IAND(s2,-8), 4), b)
b  = ISHFT( IEOR( ISHFT(s3,3), s3), -11)
s3 = IEOR( ISHFT( IAND(s3,-16), 17), b)
random_numb = IEOR( IEOR(s1,s2), s3) * 2.3283064365E-10_dp + 0.5_dp

RETURN
END FUNCTION taus88

END MODULE Ecuyer_random

MODULE mt19937
! A Fortran-program for MT19937: Real number version
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-11-26  Time: 17:09:23
! Latest revision - 5 February 2002
! A new seed initialization routine has been added based upon the new
! C version dated 26 January 2002.
! This version assumes that integer overflows do NOT cause crashes.
! This version is compatible with Lahey's ELF90 compiler,
! and should be compatible with most full Fortran 90 or 95 compilers.
! Notice the strange way in which umask is specified for ELF90.
 
!   genrand() generates one pseudorandom real number (double) which is
! uniformly distributed on [0,1]-interval, for each call.
! sgenrand(seed) set initial values to the working area of 624 words.
! Before genrand(), sgenrand(seed) must be called once.  (seed is any 32-bit
! integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.

! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Library General Public License as published by
! the Free Software Foundation; either version 2 of the License, or (at your
! option) any later version.   This library is distributed in the hope that
! it will be useful, but WITHOUT ANY WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General Public License
! along with this library; if not, write to the Free Foundation, Inc.,
! 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.

!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.

!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed

! This program uses the following standard intrinsics.
!   ishft(i,n): If n > 0, shifts bits in i by n positions to left.
!               If n < 0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.

!***********************************************************************

IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

! Period parameters
INTEGER, PARAMETER :: n = 624, n1 = n+1, m = 397, mata = -1727483681
!                                    constant vector a
INTEGER, PARAMETER :: umask = -2147483647 - 1
!                                    most significant w-r bits
INTEGER, PARAMETER :: lmask =  2147483647
!                                    least significant r bits
! Tempering parameters
INTEGER, PARAMETER :: tmaskb= -1658038656, tmaskc= -272236544

!                     the array for the state vector
INTEGER, SAVE      :: mt(0:n-1), mti = n1
!                     mti==N+1 means mt[N] is not initialized

PRIVATE
PUBLIC :: dp, sgrnd, grnd, init_genrand

CONTAINS


SUBROUTINE sgrnd(seed)
! This is the original version of the seeding routine.
! It was replaced in the Japanese version in C on 26 January 2002
! It is recommended that routine init_genrand is used instead.

INTEGER, INTENT(IN)   :: seed

!    setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
!    [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]

mt(0)= IAND(seed, -1)
DO  mti=1,n-1
  mt(mti) = IAND(69069 * mt(mti-1), -1)
END DO

RETURN
END SUBROUTINE sgrnd
!***********************************************************************

SUBROUTINE init_genrand(seed)
! This initialization is based upon the multiplier given on p.106 of the
! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.

! This version assumes that integer overflow does NOT cause a crash.

INTEGER, INTENT(IN)  :: seed

INTEGER  :: latest

mt(0) = seed
latest = seed
DO mti = 1, n-1
  latest = IEOR( latest, ISHFT( latest, -30 ) )
  latest = latest * 1812433253 + mti
  mt(mti) = latest
END DO



RETURN
END SUBROUTINE init_genrand
!***********************************************************************

FUNCTION grnd() RESULT(fn_val)
REAL (dp) :: fn_val

INTEGER, SAVE :: mag01(0:1) = (/ 0, mata /)
!                        mag01(x) = x * MATA for x=0,1
INTEGER       :: kk, y

! These statement functions have been replaced with separate functions
! tshftu(y) = ISHFT(y,-11)
! tshfts(y) = ISHFT(y,7)
! tshftt(y) = ISHFT(y,15)
! tshftl(y) = ISHFT(y,-18)

IF(mti >= n) THEN
!                       generate N words at one time
  IF(mti == n+1) THEN
!                            if sgrnd() has not been called,
    CALL sgrnd(4357)
!                              a default initial seed is used
  END IF
  
  DO  kk = 0, n-m-1
    y = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
    mt(kk) = IEOR(IEOR(mt(kk+m), ISHFT(y,-1)),mag01(IAND(y,1)))
  END DO
  DO  kk = n-m, n-2
    y = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
    mt(kk) = IEOR(IEOR(mt(kk+(m-n)), ISHFT(y,-1)),mag01(IAND(y,1)))
  END DO
  y = IOR(IAND(mt(n-1),umask), IAND(mt(0),lmask))
  mt(n-1) = IEOR(IEOR(mt(m-1), ISHFT(y,-1)),mag01(IAND(y,1)))
  mti = 0
END IF

y = mt(mti)
mti = mti + 1
y = IEOR(y, tshftu(y))
y = IEOR(y, IAND(tshfts(y),tmaskb))
y = IEOR(y, IAND(tshftt(y),tmaskc))
y = IEOR(y, tshftl(y))

! old code AL
!IF(y < 0) THEN
!  fn_val = (DBLE(y) + 2.0D0**32) / (2.0D0**32 - 1.0D0) 
!ELSE
!  fn_val = DBLE(y) / (2.0D0**32 - 1.0D0)
!END IF

! to make it (0-1) AL
IF(y < 0) THEN
  fn_val = (DBLE(y) + 2.0D0**32) 
ELSE
  fn_val = DBLE(y) 
END IF
fn_val=(fn_val+0.5d0)/4294967296.d0


RETURN
END FUNCTION grnd


FUNCTION tshftu(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,-11)
RETURN
END FUNCTION tshftu


FUNCTION tshfts(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,7)
RETURN
END FUNCTION tshfts


FUNCTION tshftt(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,15)
RETURN
END FUNCTION tshftt


FUNCTION tshftl(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,-18)
RETURN
END FUNCTION tshftl

END MODULE mt19937





module matrix_algebra
use kinds
double precision:: tol=1d-20

contains

! ----------------------------------------------------------------------


!Sacada de los programas de misztal EV. He quitado rank y he puesto n y he puesto entrada.
!    En principio, esta te devuelve la triangular inferior. Lo hemos cambiado para que 
!    te de la triangular superior y el programa sea coherente.

    SUBROUTINE CHOLS4(entrada,x)
    ! cholesky decomposition
    integer :: n,i,j,k,rank_o
    integer::rank
    double precision::diagsq
    double precision::x(:,:), entrada(:,:)

    rank_o=0
    x=entrada
    n=size(x,1)
    do i=1,n
       diagsq=x(i,i)-dot_product(x(i,1:i-1),x(i,1:i-1))
       if (abs(diagsq).lt.tol) then
         x(i,:)=0;x(:,i)=0       !zero row and column
       elseif (diagsq.lt.0) then
         print*,' Matrix not semipositive-definite, row ',i
         stop
       else
         rank_o=rank_o+1
         x(i,i)=sqrt(diagsq)
         do j=i+1,n     
           x(j,i)=(x(j,i)-dot_product(x(j,1:i-1),x(i,1:i-1)))/x(i,i)
           x(i,j)=x(j,i)
         enddo
       end if
    enddo

      ! zero upper-diagonals
      do i=1,n
         x(i,i+1:n)=0
      enddo   
    
!    Cambiamos a triangular superior 

    x=transpose(x)
      
    END SUBROUTINE CHOLS4


! ----------------------------------------------------------
      subroutine inver2(var,n)
        double precision VAR(:,:),VER(size(var,1)),VIR(size(var,1),size(var,1)),A,B
!        CALL TIME(LSEG)
!       ----------------------------------------------------------
!       LA MATRIZ CARGADA DE DECIMALES SE UTILIZA PARA DETECTAR 
!       ERRORES DE REDONDEO
!       ----------------------------------------------------------
       if(n>size(var,1)) then
         print *,'n>size of var in inver2',n,size(var,1)
         stop
       endif
       

       DO 4 I=1,N
!       INVERSION DE LA ESQUINA
         DO J=1,I-1
           VER(J)=0.0
           DO  K=1,I-1
             VER(J)=VER(J)+VAR(I,K)*VAR(K,J)
           enddo
           VAR(I,I)=VAR(I,I)-VER(J)*VAR(I,J)
         enddo
         VAR(I,I)=1.0/VAR(I,I)

!        INVERSION DEL MARGEN
         DO  J=1,I-1
           VAR(I,J)=-VAR(I,I)*VER(J)
         enddo

!        INVERSION DE LA MATRIZ BLOQUE
         DO  J=1,I-1
           DO  K=1,I-1
             VAR(J,K)=VAR(J,K)-VER(J)*VAR(I,K)
           enddo
	 enddo

!        INVERSION DEL OTRO MARGEN
         DO  J=1,N-1
           VAR(J,I)=VAR(I,J)
         enddo
4      CONTINUE

       END subroutine
    

!Sacada de los programas de misztal EV
     SUBROUTINE INVCS4(l)
! calculates inverse of LL', where L is cholesky decomposition 
      integer  ::n,i,j,rank_o
      double precision :: l(:,:)
      double precision ::w(size(l,1))
    double precision ::kkita(size(l,1),size(l,1))
    
    n=size(l,1)
    call chols4(l,kkita)     ! cholesky factorization
!      Modificado por EV
    l=kkita

    do i=1,n
      w(i:n)=0
      if (abs(l(i,i)).gt.tol) w(i)=1/l(i,i)
  ! forward substitution
      do j=i+1,n
      if (abs(l(j,j)).gt.tol) &
      w(j)=-dot_product(l(j,i:j-1),w(i:j-1))/l(j,j)
      enddo
    !backward substitution
      do j=n,i,-1
      if (abs(l(j,j)).gt.tol) &
     w(j)=(w(j)-dot_product(l(j+1:n,j),w(j+1:n)))/l(j,j)
      enddo
      l(i:n,i)=w(i:n)
      l(i,i:n)=w(i:n)
      enddo   


    END SUBROUTINE INVCS4


    

! --------------------------------------
      subroutine ginv1(a,n,m,tol,irank)
! returns generalized inverse of matrix x of size n x n declared
! as m x m. tol is working zero and irank returns the rank of
! the matrix. w is a work vector of size m,
! by rohan fernando, slightly structured by i. misztal 05/05/87

      double precision:: a(m,m),w(m),re,sum,tol
      irank=n
      do 10 i=1,n
         do 20 j=1,i-1
              re=a(i,j)
              do 20 ii=i,n
20                 a(ii,i)=a(ii,i)-re*a(ii,j)
         if (a(i,i).lt.tol) then
              a(i,i)=0.0
              do 45 ii=i+1,n
45                 a(ii,i)=0.
           irank=irank-1
           else
              a(i,i)=sqrt(a(i,i))
              do 40 ii=i+1,n
40                a(ii,i)=a(ii,i)/a(i,i)
         endif
10    continue
 
      do 100 i=1,n
         if (a(i,i).eq.0.) then
              do 150 ii=i+1,n
150                a(ii,i)=0.
           else
              a(i,i)=1.0/ a(i,i)
              do 200 ii=i+1,n
200               w(ii)=0.0
              do 300 ii=i+1,n
                  iim1=ii-1
                  re=a(iim1,i)
                  do 400 iii=ii,n
400                   w(iii)=w(iii)-a(iii,iim1)*re
                  if (a(ii,ii).eq.0.) then
                      a(ii,i)=0.
                    else
                      a(ii,i)=w(ii)/a(ii,ii)
                  endif
300           continue
          endif
100     continue
 
      do 110 j=1,n
         do 110 i=j,n
              sum=0
              do 130 ii=i,n
130                sum=sum+a(ii,j)*a(ii,i)
110           a(i,j)=sum
      do 600 i=1,n
          do 600 j=i,n
600           a(i,j)=a(j,i)

      end subroutine


end module matrix_algebra


module random_generator
  use kinds
  use matrix_algebra
  use Ecuyer_random
  use mt19937
  
  contains

!     --------------------------------------------
      subroutine unif_old(x1,u)
!      generacion de un numero uniforme u[0,1]
!      x1 es la semilla
! AL  dio problemas de ciclicidad
      implicit double precision(a-h,o-z)
      divis=2.**31.-1.
      trans=7**5
      divid=trans*x1
      lsol=int(divid/divis)
      x1=divid-lsol*divis
      u=x1/divis
      end subroutine

!     ---------------------
      subroutine unif(ix,u)
!      use kinds

!----------------------------------------------
! SUBRUTINA PARA MUESTREO DE UNIFORME- NUEVA 
! La semilla debe ser un entero comprendido entre:
! ix=2140000001      ! 0<ix<2147483647
!----------------------------------------------
      implicit none
      integer ix,k1
      double precision u
!      k1=ix/127773
!      ix=16807*(ix-k1*127773)-k1*2836
!      if(ix.lt.0)ix=ix+2147483647
!      u=ix*4.656612875E-10
      u=taus88()
!      u=grnd() !Mersenne-Twister
      end subroutine


      function fnormal() result(z)
!      generacion de un numero normal z -> n(0,1)
!      x1 es la semilla
      implicit double precision(a-h,o-z)
      double precision:: z,u1,u2
      call unif(0,u1)
      call unif(0,u2)
      z=((-2.*log(u1))**0.5)*cos(2.*3.1416*u2)
      end function

      subroutine normal(x1,z)
!      generacion de un numero normal z -> n(0,1)
!      x1 es la semilla
      implicit double precision(a-h,o-z)
    integer x1
      double precision:: z,u1,u2
      call unif(x1,u1)
      call unif(x1,u2)
      z=((-2.*log(u1))**0.5)*cos(2.*3.1416*u2)
      end subroutine

!     ================Nuevas==============

! -----------------------------------------------------
      subroutine gen_trunc_normal_old(cota1,cota2,z,x1)
!      generacion de un numero normal z -> normal truncada entre cota1 y cota2
!     (cotas estandarizadas a varianza=1)
!      x1 es la semilla
!     Metodo brutal a base de prueba y error 
!     Existe un metodo mas fino en Sorensen, p.146
!     AL 15/4/04
      implicit double precision(a-h,o-z)
      double precision:: z,xnor
      double precision:: cota1,cota2
      integer i,x1

      i=0
      do 
        i=i+1
!        if (i.gt.10000000) then
!          print *,i
!          stop
!        endif
!        z=xnor(x1)
        call normal(x1,z)
        if ((z.gt.cota1).and.(z.lt.cota2)) goto 101
      enddo
 101  continue
      end  subroutine

! ----------------------------------------------------------------------
      subroutine chin(se,ne,x,u)
!      se = varianza a priori.
!      ne = grados de credibilidad.
      implicit double precision (a-h,o-z)
      double precision ne
      integer x
1     call unif(x,u1)
      rg=se*14/(ne**.5)
      t=2*rg*u1+(se-rg)
      if (t.lt.0) goto 1
      o=(ne/(ne+2))*se
      chiot=-((ne+2)/2)*log(o)+(-(ne*se/2)/o)
      chitt=-((ne+2)/2)*log(t)+(-(ne*se/2)/t)
      ref=dexp(chitt-chiot)
      call unif(x,u1)
      if (u1.lt.ref) then
        u=t
      else
        goto 1
      endif

      end subroutine


      subroutine wish(nrank,se,ve,ne,x1)
       implicit double precision(a-h,o-z)
       double precision se(:,:),ve(:,:)
       double precision t(size(se,1),size(se,1)),a(size(se,1),size(se,1)),l(size(se,1),size(se,1))
       double precision b(size(se,1),size(se,1))
       integer x1
!      double precision se(nrank,nrank),ve(nrank,nrank)
!      double precision t(nrank,nrank),a(nrank,nrank),l(nrank,nrank)
!       double precision b(nrank,nrank)
         ia=ne/2

       do 1 i=1,nrank
           au=real(ne-i+1)/2.
           fio=gamdev(au,x1)
         t(i,i)=sqrt(2*fio)
         do 2 j=1,i-1
!           u=xnor(x1)
           call normal(x1,u)
            t(i,j)=0.
            t(j,i)=u
2         continue
1       continue
       do 3 i=1,nrank
          do 4 j=1,nrank
             a(i,j)=0.
             do 5 k=1,nrank
             a(i,j)=a(i,j)+t(i,k)*t(j,k)
5             continue
4          continue
3       continue
       do 6 i=1,nrank
         do 6 j=1,nrank
           l(i,j)=0.
6      continue
       call chols4(se,l)
       do 30 i=1,nrank
          do 40 j=1,nrank
             b(i,j)=0.
             do 50 k=1,nrank
               b(i,j)=b(i,j)+l(k,i)*a(k,j)
50           continue
40        continue
30     continue
       do 31 i=1,nrank
          do 41 j=1,nrank
             ve(i,j)=0.
             do 51 k=1,nrank
               ve(i,j)=ve(i,j)+b(i,k)*l(k,j)
51             continue
41        continue
31     continue
       do 101 i=1,nrank
         do 101 j=1,nrank
           ve(i,j)=ve(i,j)/ne
101    continue

       end subroutine


    


! ---------------------------------------------------      
      subroutine inv_con_wish(ncar,se,ve,ne,x1)
!      muestrea una wishart invertida para una suma de cuadrados (NO para su inversa)
!      condicionada a que un elemento 
!      (el del car\E1cter umbral) sea igual a 1.
!      Teor\EDa de Korsgaard, de 
!       http://www.csc.fi/ttn/ccb99/articles/IKorsgaard.pdf
!      Solo trabaja con el ultimo caracter
!      AL, 15/4/2004
!      10/6/04 -> corregido error en calculo de meant2
!      (no hacia multiplicacion matricial)

      double precision:: se(:,:),ve(:,:)
      integer x1
      integer ne,i,j,ncar
!     auxiliares
      double precision:: se_old(size(se,1),size(se,1))
      double precision:: V11(size(se,1),size(se,1)),V11_old(size(se,1),size(se,1))
      double precision:: t2(size(se,1)),meant2(size(se,1)),vart2(size(se,1),size(se,1))

! se como en el paper (sumas de cuadrados y no sumas de cuadrados entre el n de datos)
      do i=1,ncar
        do j=1,ncar
          se(i,j)=se(i,j)*real(ne)
        enddo
      enddo
! lo invierto para q sea como en el paper
      call inver2(se,ncar)
! y lo guardo
      do i=1,ncar
        do j=1,ncar
          se_old(i,j)=se(i,j)
        enddo
      enddo

! 1- V11
      call wish(ncar-1,se,V11,ne-1,x1)
      do i=1,ncar-1
        do j=1,ncar-1
! por la parametrizacion de la subrutina Wishart del paper donde E(W|se,ne)=se*ne
! wish da E(wish)=se
! por tanto lo reescalo
          V11(i,j)=V11(i,j)*real(ne-1)
          V11_old(i,j)=V11(i,j)
        enddo
      enddo
! inv(V11)
      call inver2(V11,ncar-1)

! 2- t2
! inv(se)
      call inver2(se,ncar)
      do i=1,ncar-1
        do j=1,ncar-1
! para sacar se_sup(22)
          vart2(i,j)=V11(i,j)*(1.d0/se(ncar,ncar))
        enddo
      enddo

! se de nuevo bien
      call inver2(se,ncar)
! para sacar inv(se_11)
      call inver2(se,ncar-1)
      do i=1,ncar-1
        meant2(i)=0.d0
        do j=1,ncar-1
          meant2(i)=meant2(i)+se(i,j)*se_old(j,ncar)
        enddo
      enddo

      call mvn(meant2,vart2,t2,ncar-1,x1)
      
      do i=1,ncar-1
        do j=1,ncar-1
          ve(i,j)=V11(i,j)+t2(i)*t2(j)
        enddo
        ve(i,ncar)=-t2(i)
        ve(ncar,i)=ve(i,ncar)
      enddo
      ve(ncar,ncar)=1.d0

      end subroutine
      
!     -------------------------------------------------------
      subroutine inv_con_wish_multiple(ncar,nres,se,ve,ne,x1)
!  muestrea una wishart invertida para una suma de cuadrados (NO para su inversa)
!  condicionada a que una submatriz 
! (la de los caracteres umbrales) sea igual a I.
! Teor\EDa de Korsgaard, de http://www.csc.fi/ttn/ccb99/articles/IKorsgaard.pdf
! AL, 10/6/04
! Ultima correccion 20/6/04
! 

! Problema2: puede producir matrices no positivas-definidas... no se si esto 
! choca con la teoria o no.
! Modified 6/7/07 to tale out "maxcar"
      ! be aware that most of this sizes anre not correct and make use of inver2(a,n) inverting n just in 
      ! its nxn upper left matrix

      double precision:: se(:,:),ve(:,:) 
! se -> sigma
      integer x1
      integer ne,i,j,ncar,k1,k2
!     auxiliares
      double precision:: se_old(size(se,1),size(se,1)),se_221(size(se,1),size(se,1))
      double precision:: V11(size(se,1),size(se,1)),V11_old(size(se,1),size(se,1))
      double precision:: t2((ncar-nres)*nres),meant2((ncar-nres)*nres)
      double precision:: t2new(size(se,1),size(se,1))
      double precision:: vart2((ncar-nres)*nres,(ncar-nres)*nres)
      integer pos1,pos2
      integer ncont,nres
! ncont -> n caracteres continuos + umbrales con var. residual estimada
!          (p.ej. con varias categorias)
! nres numero de caracteres umbral con var res restringida a 1
      ncont=ncar-nres

! se como en el paper (sumas de cuadrados y no sumas de cuadrados 
! entre el n de datos)
      do i=1,ncar
        do j=1,ncar
          se(i,j)=se(i,j)*real(ne)
        enddo
      enddo
! lo invierto para q sea como en el paper
      call inver2(se,ncar)
! y lo guardo
      do i=1,ncar
        do j=1,ncar
          se_old(i,j)=se(i,j)
        enddo
      enddo

! 1- V11
      call wish(ncont,se,V11,ne-nres,x1)
      do i=1,ncont
        do j=1,ncont
! por la parametrizacion de la subrutina Wishart del paper donde E(W|se,ne)=se*ne
! mientras que la subrutina wish da E(wish)=se
! por tanto lo reescalo
          V11(i,j)=V11(i,j)*real(ne-nres)
          V11_old(i,j)=V11(i,j)
        enddo
      enddo
! inv(V11) -> v_11exp-1
      call inver2(V11,ncont)

! 2- t2
! inv(se)-> sigma-1
      call inver2(se,ncar)
! extraigo de sigma-1 el bloque correspondiente a se_22.1    
      do i=1,nres
        do j=1,nres
          se_221(i,j)=se(ncont+i,ncont+j)
        enddo
      enddo
! y lo invierto -> se_22.1
      call inver2(se_221,nres)

      do i=1,ncont
        do j=1,ncont
          do k1=1,nres
            do k2=1,nres
              pos1=(i-1)*nres+k1
              pos2=(j-1)*nres+k2
              vart2(pos1,pos2)=V11(i,j)*se_221(k1,k2)
            enddo
          enddo
        enddo
      enddo

! se de nuevo bien ->sigma
      call inver2(se,ncar)
! para sacar inv(se_11)
      call inver2(se,ncont)
      do i=1,ncont
        do k1=1,nres
          pos1=(i-1)*nres+k1
          meant2(pos1)=0.d0
          do j=1,ncont
            meant2(pos1)=meant2(pos1)+se(i,j)*se_old(j,ncont+k1)
          enddo
        enddo
      enddo

      call mvn(meant2,vart2,t2,nres*ncont,x1)
      
!      recoloco t2
      do i=1,ncont
        do k1=1,nres
          pos1=(i-1)*nres+k1
          t2new(i,k1)=t2(pos1)
        enddo
      enddo

! v11-1+t2t2'
      do i=1,ncont
        do j=1,ncont
          ve(i,j)=V11(i,j)
          do k1=1,nres
            do k2=1,nres
              ve(i,j)=ve(i,j)+t2new(i,k1)*t2new(j,k2)
            enddo
          enddo
        enddo
      enddo

! -t2
      do i=1,ncont
        do k1=1,nres
          ve(i,ncont+k1)=-t2new(i,k1)
          ve(ncont+k1,i)=ve(i,ncont+k1)
        enddo
      enddo

! Im4
      do i=1,nres
        do j=1,nres
          if (i.eq.j) then
            ve(ncont+i,ncont+j)=1.d0
          else 
            ve(ncont+i,ncont+j)=0.d0
          endif
        enddo
      enddo
      
! Chequeo si la inversa no tiene negativos
      do i=1,ncar
        do j=1,ncar
          se_old(i,j)=ve(i,j)
        enddo
      enddo
      call inver2(se_old,ncar)
      do i=1,ncar
          if (se_old(i,i).lt.0.d0) then
            print*,'negativo en inv(ve)',i,i,se_old(i,i)
            do j=1,ncar
              print *,(ve(j,k),k=1,ncar)
            enddo
          endif
      enddo

      end subroutine


! -------------------------------------------------
      subroutine inv_wish_special(rank,se,ve,ne,x1)
! Devuelve wishart invertidas para matrices con elementos operacionalmente 0
! Para no tener problemas numericos
! AL
      double precision:: se(:,:),ve(:,:),zero
      integer ne,rank,kk,x1
! Cero operacional
      zero=1e-15     

      call ginv1(se,rank,size(se,1),zero,kk)
      call wish(rank,se,ve,ne,x1)
      call ginv1(ve,rank,size(se,1),zero,kk)

      end subroutine


! ------------------------------
      subroutine mvn(a,v,b,n,x1)
! Genera normales multivariadas b~MVN(a,v) de rango n con semilla x1
! A partir de una subrutina de Luis Varona
      integer i,n,j,x1
      double precision a(:),v(:,:),cv(size(v,1),size(v,1)),b(n)
      double precision u(n),xnor
      call chols4(v,cv) 
      do i=1,n
        call normal(x1,u(i))
!        u(i)=xnor(x1)
      enddo
      do i=1,n
        b(i)=a(i)
        do j=1,i
          b(i)=b(i)+u(j)*cv(j,i)
        enddo
      enddo
      end subroutine


!---- Subrutinas para obtener normales truncadas; principalmente de programas
!       de Misztal

! ---------------------------------
      function normal_invcdf(p) 
      ! return inverse of CDF, i.e., such x: p=cdf(x)
      !   resticted to |invcdf(x)|<10

!      use kinds
      double precision:: p,x,eps,normal_invcdf
      integer i

      if (p.lt.0. .or. p .gt. 1.) then
         print*,'normal_invcdf: arguments outside of 0..1'
         stop
      endif

      eps=10

      x=0
      do i=1,50
         if (normalcdf(x) .lt. p) then     
           x=x+eps
         else
           x=x-eps
         endif    
         if (abs(x) .gt. 10.1) exit
         eps=eps/2.
      enddo      
      normal_invcdf=x
      end function

! ---------------------------
      function normalcdf(x)
      ! returns cumulative normal density function
      implicit none
      double precision:: normalcdf,x,alnorm,q,pdf
      !
      !normalcdf=alnorm(x,.false.)
       call  NORMP(x, normalcdf, Q, PDF)

      end function


!-----------------------------------
      SUBROUTINE NORMP(Z, P, Q, PDF)
!
!    Normal distribution probabilities accurate to 1.e-15.
!    Z = no. of standard deviations from the mean.
!    P, Q = probabilities to the left & right of Z.   P + Q = 1.
!       PDF = the probability density.
!
!       Based upon algorithm 5666 for the error function, from:
!       Hart, J.F. et al, 'Computer Approximations', Wiley 1968
!
!       Programmer: Alan Miller
!
!    Latest revision - 30 March 1986
!
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DATA P0, P1, P2, P3, P4, P5, P6/220.2068679123761D0, &
           221.2135961699311D0, 112.0792914978709D0, &
           33.91286607838300D0, 6.373962203531650D0, &
          .7003830644436881D0, .3526249659989109D-01/, &
           Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7/440.4137358247522D0, &
           793.8265125199484D0, 637.3336333788311D0, &
           296.5642487796737D0, 86.78073220294608D0, &
           16.06417757920695D0, 1.755667163182642D0, &
           .8838834764831844D-1/, &
           CUTOFF/7.071D0/, ROOT2PI/2.506628274631001D0/
!
      ZABS = ABS(Z)
!
!      |Z| > 37.
!
      IF (ZABS .GT. 37.D0) THEN
        PDF = 0.D0
        IF (Z .GT. 0.D0) THEN
          P = 1.D0
          Q = 0.D0
        ELSE
          P = 0.D0
          Q = 1.D0
        END IF
        RETURN
      END IF
!
!      |Z| <= 37.
!
      EXPNTL = EXP(-0.5D0*ZABS**2)
      PDF = EXPNTL/ROOT2PI
!
!      |Z| < CUTOFF = 10/sqrt(2).
!
      IF (ZABS .LT. CUTOFF) THEN
        P = EXPNTL*((((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS + &
             P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS + &
             Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS + &
             Q0)
!
!      |Z| >= CUTOFF.
!
      ELSE
        P = PDF/(ZABS + 1.D0/(ZABS + 2.D0/(ZABS + 3.D0/(ZABS + 4.D0/ &
         (ZABS + 0.65D0)))))
      END IF
!
      IF (Z .LT. 0.D0) THEN
        Q = 1.D0 - P
      ELSE
        Q = P
        P = 1.D0 - Q
      END IF
      END subroutine

! ---------------------------
      function dnormal(x)
      !
      ! Returns probability density of normal density function
      !
      implicit none
      double precision:: dnormal,x
      double precision:: pi
      
      pi=3.14159265358979
      dnormal=1.d0/sqrt(2.d0*pi)*dexp(-x*x/2.d0)
      end function

      subroutine gen_trunc_normal(a,b,x,x1)
! Generates normal distribution N(0,1) truncated to <a,b>

      double precision:: a,b,cdfa,cdfb
      double precision:: un,x
      integer x1


      if (b .lt. a) then
         Print*,'GEN_TRUNC_NORMAL: bound b =',b,' < boound a =',a
         stop
      endif   

      cdfa=normalcdf(a)
      cdfb=normalcdf(b)
1     continue        
      call unif(x1,un)  !uniform random number generator UN(0,1)
      x=normal_invcdf(cdfa+(cdfb-cdfa)*un)       
      if ((x.lt.a).or.(x.gt.b)) then
        print *,a,x,b
	print *,'I can''t sample this truncated normal, subroutine gen_trunc_normal'
!        goto 1
      endif
      end subroutine
      
! ----------------------------------
       FUNCTION XNORMD(PROB)
!       SOLO DEFINIDA PARA PROB EN EL ESPACIO 0-.5
       IMPLICIT double precision (A-H,O-Z)
       double precision,save:: limit=tiny(1d0)
       if(prob<limit) prob=limit !AL to avoid log(0)
       T=SQRT(LOG(1./(PROB*PROB)))
       XNORMD=T-(2.515517+T*(.802853+T*.010328))/(1.+T* &
       (1.432788+T*(.189269+T*.001308)))
       END function

! -------------------------------------

      function xnor_i(x1) result(u)
      integer :: x1
      double precision:: u
      
      u=fnormal()
      end function
      
      
      FUNCTION XNOR(x1)

       IMPLICIT double precision (A-H,O-Z)
       integer x1
       double precision:: xnor

!       CALCULA LA ABCISA (X) PARA CUALQUIER INCIDENCIA (P) 
! Genera numeros de distribucion normal mediante transformacion a partir
! de una uniforme (x1 es el valor de esa uniforme)
! Segun Luis no es tan fina como (normal) pero es mejor en el modelo umbral
! para muestrear normales truncadas "a pelo"
! Modificada AL para usar unif() ya que la otra subrutina dio problemas
! de ciclicidad  

! Generacion de muestra de la uniforme

!       divis=2.**31.-1.
!       trans=7.**5.
!       divid=trans*x1
!       lsol=int(divid/divis)
!       x1=divid-lsol*divis
!       prob=x1/divis

       call unif(x1,prob)

!      Transformacion
       IF(PROB.LT..5) THEN
          XNOR=-XNORMD(PROB)
       ELSE IF(PROB.GT..5)THEN
          XNOR=XNORMD(1.-PROB)
       ELSE
          XNOR=0.
       ENDIF
       END function

!     ----------------------------------------------------------------
      function gamdev(a,x1)
       implicit double precision(a-h,o-z)
       integer x1
       if(a.lt.1)then
         print *,'illegal alpha parameter in gamdev'
       endif
       if(a.lt.6)then
         x=1.
           ia=int(a)
         do 11 j=1,ia
            call unif(x1,u)
            x=x*u
11        continue
         x=-log(x)
       else
1         call unif(x1,u)
         v1=2.*u-1.
         call unif(x1,u)
         v2=2.*u-1.
         if(v1**2+v2**2.gt.1.)goto 1
         y=v2/v1
           am=a-1
         s=sqrt(2.*am+1.)
         x=s*y+am
         if (x.le.0.)goto 1
            e=(1.+y**2)*exp(am*log(x/am)-s*y)
         call unif(x1,u)
         if (u.gt.e)goto 1
       endif
       gamdev=x
       end function



end module random_generator



module aux1

contains

      subroutine printtime
      INTEGER  DATE_TIME (8)
      CHARACTER(LEN = 12) REAL_CLOCK (3)

      CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), &
                      REAL_CLOCK (3), DATE_TIME)

      print *,real_clock(1)(7:8),'/',real_clock(1)(5:6),&
      '/',real_clock(1)(1:4),' ', &
      real_clock(2)(1:2),':',real_clock(2)(3:4),':',real_clock(2)(5:6)

      end subroutine

      subroutine crearotulo(ncar,nefani,nper,rotulo,pos1,tipomodelo)

      dimension rotulo(100)
      integer ncar,nper,i,j,k,pos1
      character*100 rotulo
      character*4 auxrotulo
      character*2 auxrotulo2 
      character(40) :: tipomodelo    

      POS1=0
      do i=1,ncar*nefani
        do j=i,ncar*nefani
          pos1=pos1+1
          write(auxrotulo,'(2i2.2)')i,j
          rotulo(pos1)=('vara_'//auxrotulo)
        enddo
      enddo
      do k=1,nper
        do i=1,ncar
          do j=i,ncar
            pos1=pos1+1
            write(auxrotulo,'(2i2.2)')i,j
            write(auxrotulo2,'(i2.2)')k
            rotulo(pos1)=('varp'//auxrotulo2//'_'//auxrotulo)
          enddo
        enddo
      enddo
      do i=1,ncar
        do j=i,ncar
          pos1=pos1+1
          write(auxrotulo,'(2i2.2)')i,j
          rotulo(pos1)=('vare_'//auxrotulo)
        enddo
      enddo
      if ((tipomodelo=='animal').and.(nefani==1)) then
	do i=1,ncar
	  pos1=pos1+1
	  write(auxrotulo,'(i2.2)')i
	  rotulo(pos1)='h2_'//auxrotulo
	enddo
      endif
      return
      end subroutine

      function sqrt1(a)
!     returns square root of a if a>0, otherwise 0
      double precision:: a,sqrt1
      if(a.le.0.d0) sqrt1=0.d0
      if(a.gt.0.d0) sqrt1=sqrt(a)
      end function

end module aux1

module aux_options
! from IM's blupps1.f90
implicit none
character(len=20):: xc(3)
integer:: nitem
real:: xcnum(3)

contains

  subroutine nums2(a,n,x,xc)
  ! separates array a into items delimited by blanks. character elements are
  ! put into optional character vector xc, decoded numeric values
  ! into optional real vector x, and n contains the number of items. The
  ! dimension of x and xc can be lower than n.
  ! A modification of nums() from f77 to f90
  ! Now accepts real numbers
  ! 2/23/2000

  character (*)::a
  character (*),optional::xc(:)
  real,optional::x(:)
  integer::n,curr,first,last,lena,stat,i

  curr=1;lena=len(a);n=0

  do
    ! search for first nonspace
    first=0
    do i=curr,lena
       if (a(i:i) /= ' ') then
          first=i
          exit
       endif
    enddo
    if (first == 0) exit


    ! search for first space
    curr=first+1
    last=0
    do i=curr,lena
        if (a(i:i) == ' ') then
          last=i
          exit
        endif
    enddo

    if (last == 0) last=lena

    n=n+1
    if (present(xc)) then
       if (size(xc) >= n) then
             xc(n)=a(first:last)
           else
             stop "size of xc in nums2 too small"
       endif
    endif
    if (present(x)) then
       if (size(x) >= n) then
              read(a(first:last),'(f12.0)',iostat=stat)x(n)
              if (stat /=0) x(n)=0
           else
              stop "size of x in nums2 too small"
       endif
    endif

    curr=last+1
  enddo
  end subroutine
  subroutine getoption_unit(io_unit,name,n,x,xc)
  ! In unit io_unit, locates line:
  ! OPTION name str1 str2
  ! where str1, str2 are strings separated by spaces.
  ! Then, it assigns: xc(1)=str1,  xc(2)=str2,...
  ! and attempts to decode strings into real values: x1=value(str1),....
  !
  ! n contains the number of strings.  x and xc are optional and their
  ! dimensions may be smaller than n in which case some strings/values are
  ! not stored.
  !
  ! Upon exit, unit io_p=40 points to line next to the one located.
  !
  ! If the line cannot be located, n=-1
  !

  character (*)::name
  integer::n,io_unit
  real,optional::x(:)
  character (*),optional::xc(:)

  real::x1(100)
  integer::stat,m
  character (50)::xc1(100)
  character (200)::a

  rewind io_unit
  n=-1

  do
     read (io_unit,'(a)',iostat=stat)a
     if (stat /= 0) exit
        call nums2(a,m,x1,xc1)
        if ( m>1 .and. xc1(1) == 'OPTION' .and. xc1(2) == name) then
           n=m-2
           if (present(xc)) xc=xc1(3:size(xc)+2)
           if (present(x)) x=x1(3:size(x)+2)
           exit
        endif
  enddo
  end subroutine
  
end module




!        ----------
         Program TM      !Threshold model
!        ----------

! 19/4/04
!    Programa para analizar un caracter categorico y varios caracteres lineales
!    permitiendo modelos diferentes para cada caracter, multiples efectos permanentes
!    y multiples efectos animales correlacionados (ej. efectos maternos). 
! 
!    Basado en un original multicaracter de Luis Varona, 
!    comentado y modificado por Andres Legarra y Evangelina Lopez de Maturana

!  15/05/04 A\F1ado multiples categorias (fenotipo=1 a n, cambiable)
!  15/05/04 Cambiado el generador de numeros aleatorios uniformes para evitar problemas de ciclicidad
!  10/06/04 Implementado algoritmo de Korsgaard para multiples binarios
!           Implementado metodo mas elegante de pelear con binarios o 
!           de varias categorias. 
!           Varias -> umbral1=0, umbral2=1
!           binarios -> umbral1=0, vare=1
!  19/06/04 Generalizado a varios categoricos 
!           Orden: lineales, policotomicos, dicotomicos
!  20/06/04 Chequeado algoritmo Korsgaard de muestreo de var(e)
!           Aumento de datos para faltantes en umbrales
!           ~N(Xb+Zu,R) -> como generamos la liability,
!           no implica a los umbrales por tanto es EXACTAMENTE
!           IGUAL a los caracteres lineales 
!  05/07/04 Inicializados sum y sum de cuadrados de componentes de varianza
!           (sintaxis f90)
!  20/07/04 Correcccion en sintaxis de write (no compilaBA CON nag) 
!           Liability para missing values limitada entre -999 y +999
!  22/07/04 Comienzo modificacion para multiples efectos animales
!           (directo + materno, etc). Solo considera 1 poblacion
!           (o sea un A-1). Si son varias poblaciones se deberia crear un mismo 
!           fichero de pedigri conjunto... con 2 grupos geneticos por cada uno 
!           e ignorar la correlacion genetica (o hacerla 0) ya que no estimable
!           Hay que meter una opcion de que la correlacion sea 0?
!           Ademas meto ef. animal variable
!           Almacenamiento de las ecuaciones en zhz:
!           X'X  X'Z1  X'Z2
!           Z1'X A-1   Z1'Z2
!           Z2'X Z2'Z1 A-1
!  26/08/04 Subrutina printtime para imprimir hora (f90)
!           Imprime cov ef permanente y rotulo
!  10/10/04 Evangelina LdM introduce modelo macho
!           Para el modelo macho hay que meter una palabra clave adicional
!           'macho'. Si es modelo macho NO se mete grupo genetico (se deja 0
!           en animales desconocidos)
!  21/10/04 Corregido bug en generacion de liability (e~MVN(e,R) mal colocado)
!  22/11/04 Corregido error en poner a 0 grupos geneticos (lo hacia cuando no los habia)
!  29/11/04 EV modifica las subrutinas chols4 (de Misztal), etc para trabajar con matrices de
!           covarianzas semipositivas definidas. Ej: permanente o efectos maternos 
!           para solo un caracter. Existen 2 subrutinas de inversion: invcs4 (de programas
!            de Misztal y ginv1 (id, de Rohan Fernando). Se utiliza esta ultima.
!  07/12/04 convertido a f90 (estilo, allocate, uso de modulos y poco mas; hay que hacer
!           bien los "allocate", cambiar las entradas de las subrutinas, etc)
!  15/12/04 Declaradas e inicializadas todas las variables y todo es real(r8).
!  16/12/04 Existian problemas de ciclicidad en problemas complejos (da+fp)
!           ya que si tol=1d-10 el programa ciclaba, y si tol < que eso
!           existian problemas numericos (vara=0, NaN, etc). Hemos cambiado el sistema
!           de muestreo segun modelo a b|resto =0.d0 si el efecto NO esta en el modelo
!           para ese caracter (antes era b|resto=N(1E-10,1E-10) ). 
!  20/12/04 Sigue habiendo problemas de ciclicidad. Incluido un generador de numeros 
!           aleatorios de L'Ecuyer sacado de la web de Alan Miller 
!  20/12/04 El programa est\E1 preparado para continuar en el caso de que el programa finalice
!           antes de lo previsto(corte de luz) o en el caso de que se quiera aumentar el n\FAmero
!           de iteraciones una vez acabado un proceso. 
!  10/01/05 Incluida subrutina pos_def de I. Misztal para chequear que las matrices de covarianzas 
!           sean semipositivas definidas. Si no lo son, hace un "bending": convierte el menor 
!           eigenvalue en 1d-6 (esto choca con el espiritu del Gibbs sampling, pero en la practica
!           no se observan problemas)
!           Esta subrutina necesita el modulo denseop q a su vez necesita kinds y lapack90
!           que van en ficheros aparte
!           AL: la verdad es q el problema solo existe cuando los analisis son muy muy problematicos
!  18/02/05 Ahora el programa, adem\E1s de permitir missing values, tambi\E9n admite datos censurados.
!           Estos se codifican con numeros negativos y se muestran de una normal truncada a partir
!           del punto de censura
!  31/03/05 Reading of initial variances in the parameter file
!           Reading of number of iterations, burn-in, and thin interval from par file
!           Choice of Task: Variance component estimation (VCE) or breeding value prediction
!           with known components of variance (BLUP)
!  11/04/05 Corrected small bug reading residual variances
!  14/04/05 real(8) is non portable syntax neither is real*8. 
!           Options: use DOUBLE PRECISION instead, or declare USE KINDS and real(r8)
!           kk and kk2 changed to logical (no speed improvement, though :( )
!  04/11/05 moved to INRA-SAGA
!           some changes in order to compile under NAG f95: real(r8) -> double precision or real(r8) and use kinds
!  06/04/06 Output in english; desactivated continuation option; manual written; optional limit for liability in binary traits
!  14/04/06 handles censored categorical traits Beware: It has no sense when the lower bound is the 1st category (same as missing value!)
!           So useful e.g. for 3 categories and values known to be >=2
!  22/05/06 Bug in writing solutions.txt (found by Jorge Urioste)
!  23/05/06 It did not invert variance components after each iteration of BLUP (e.g., took inverses for real parameters)
!           Thanks to JU again
!  11/09/06 Messages in english
!           Do not print h2 and rg for sire or multiple animal effects model
!  06/07/07 Rewrote using reorganized modules
!  07/07/07 All storage is dynamic (no recompilation needed!)


!-- Por hacer:
!     test de cadenas paralelas
!     la continuacion no esta bien hecha ya que en lugar de leer s1 s2 s3 lee la semilla x1
!     (que ya no se usa)
!     X'R-1y se arma para todos los efectos independientemente de si estan en el modelo o no
!     no es problematico pq la ecuacion no se usa, pero hay que corregirlo para evitar problemas
!     si se hiciera muestreo por bloques,etc
!     covariables
!     Foulley advice is to use two different RNG in different executions
!     
! $Id: tm.f90,v 1.7 2007/07/06 07:45:36 legarra Exp $
! $Author: legarra $
! $Date: 2007/07/06 07:45:36 $
! $Name:  $
! $Log: tm.f90,v $
! Revision 1.7  2007/07/06 07:45:36  legarra
! mistake
!
! Revision 1.6  2007/07/06 07:36:43  legarra
! Elimination of maxfac, etc
!
! Revision 1.5  2007/07/06 06:43:40  legarra
! Finished cleaning of matrix algebra and random modules
!
! Revision 1.4  2007/07/05 13:10:40  legarra
! First rewrite of TCM as to put everything in modules
!
! Revision 1.3  2007/07/05 11:31:05  legarra
! just to check $Id: tm.f90,v 1.7 2007/07/06 07:45:36 legarra Exp $ and related stuff
!
! $RCSfile: tm.f90,v $
! $Revision: 1.7 $
! $Source: /home/legarra/cvs_mine/TCM/tm.f90,v $


! 10/03/08 Bug corrected in missing record for categorical trait
! 11/03/08 Implemented Mersenne-Twister RNG (see subroutine unif()) but not used as it is slower
!          Note these times:
!          xnor, Ecuyer -> 4.697 s for 10^8 normal samples
!          xnor, Mersenne-Twister -> 6.516
!          fnormal, Ecuyer -> 10.255
!          fnormal, Mersenne-Twister -> 15.610
! I do not know if this applies to the program as most of the time is spent looking up the linked list up and down
! 17/10/08 Merged with Evangelina version which considers covariables.
! 02/08/11 Optional random seeds


! use denseop temporary disabled AL
use kinds
use modelo
use listaligada
use matrix_algebra
use random_generator
use aux1
use aux_options
use vers

implicit none

!
real(r8):: uu, ee, landa
!
real(r8),allocatable::  b(:,:)
real(r8):: u,var
real(r8),allocatable:: vara(:,:),vare(:,:),varp(:,:,:),vark(:)
real(r8),allocatable:: se(:,:),ve(:,:),va(:,:),sa(:,:)
real(r8):: ne
real(r8),allocatable:: y(:,:)
real(r8),allocatable:: reg(:,:)
real(r8),allocatable:: thr(:,:)
real(r8):: cota1,cota2
real(r8):: total
real(r8),allocatable:: mint(:,:),maxt(:,:)
real(r8):: xm


!ev_mad
real(r8),allocatable::valor(:,:)
!
integer:: npdve,npdvp,npdva ! numero de veces q Ve Vp y Va son no semipositivo definido
logical::  stat
integer::  nfac
!ev_mad
integer:: ncov
!
integer::  ngrup
integer::  nefani
integer::  i,j,k,ijk,kji,jk,ji,ki
logical::  kk
integer::  nrow
integer::  mm
integer::  nanim
integer::  nangru
integer::  mami
integer::  ia,ip,im
integer:: imue,lap,icad
integer::  is,iss,imgs
integer::  ndat1
integer::  iplace
integer::  ilev
integer::  irango
integer::  l
integer :: na
integer::  itercor, itermas, respuesta
integer::  ncar,nsegm,nyear
integer::  x1
!     AL 
integer::  presenthr,nper,iper
integer::  ncarumb,ncont,poscar,nres,maxthr
integer::  pos1,pos2,pos3,iefani,iefani2
integer::  posani
integer::  ia1,ip1,im1,IS1,ISS1,IMGS1,ig,iy
character *2 aux,aux2
logical,allocatable :: kk2(:,:)
integer,allocatable :: ifac(:),nfix(:,:),mfac(:),ind(:,:),nthr(:), &
                       posvara1(:,:),posvara2(:,:),iped(:,:),ny(:),ndesc(:),nwei(:),nsel(:)
real(r8),allocatable :: nzz(:),dia(:),yy(:,:),sumy(:,:,:,:),sumwei(:,:,:),sumsel(:,:,:)

real(r8),allocatable:: sumvara(:,:),ssvara(:,:),sumvarp(:,:,:),ssvarp(:,:,:), &
                     sumvare(:,:),ssvare(:,:),bv(:,:,:,:)
real(r8),allocatable:: ssumyeee(:,:,:,:),sssumyeee(:,:,:,:)

real(r8),allocatable:: cora(:,:),corp(:,:,:),core(:,:),h2(:)
real(r8),allocatable:: sumcora(:,:),sumcore(:,:),sumcorp(:,:,:),sscora(:,:), &
                     sscore(:,:),sscorp(:,:,:)
real(r8),allocatable:: sumb(:,:),ssb(:,:),xmend(:,:,:)
real(r8),allocatable:: varyear(:,:,:),svaryear(:,:,:),ssvaryear(:,:,:),ssvar(:,:,:)
real(r8),allocatable:: differ(:,:,:),intens(:,:,:)
real(r8),allocatable:: sdiffer(:,:,:),ssdiffer(:,:,:)
real(r8),allocatable:: sintens(:,:,:),ssintens(:,:,:)
real(r8):: x,x2,xmean,sd
real:: seeds(3)=0d0
character(len=40) fichero,datos,genea
character(len=100) rotulo(100)
character(len=40) tipomodelo,recursividad
character(10)::kkk,task
logical :: VCE=.true.
!EV
real(r8),allocatable:: d(:),v(:,:) !Eigenvalues y eigenvectors
real(r8) :: liabilitybound=4.d0 !other option is 4.d0 to constrain it

! to get contrasts
integer:: efftest,sizecontrast,io
integer,allocatable:: levefftest(:)
character(5):: contrast


call print_version(start=.true.)


npdve=0.d0; npdvp=0.d0; npdva=0.d0
itercor=0; itermas=0;respuesta=0
total=0.d0
u=0.d0
presenthr=0
iper=0
poscar=0
pos1=0
pos2=0
iefani=0
iefani2=0
posani=0
ia1=0
ip1=0
im1=0
is1=0
iss1=0
imgs1=0
nfac=0;ngrup=0;nefani=0;l=0;i=0;j=0;ijk=0;kji=0;jk=0;ji=0;mm=0;nanim=0;nangru=0
ia=0;ip=0;im=0
is=0;iss=0;imgs=0
iplace=0
irango=0
ilev=0
na=0
! Pone a 0 
NPLACE=0
!EV



! --> Ficheros de salida
open(15,file='samples.txt',position='append')
open(16,file='thresholds.txt',position='append')
open(20, file='samplesFE.txt',position='append')


!-----------------------------------------
! --> Parameter file reading <--
!-----------------------------------------

   
print *,'Parameter file?'
fichero='param.txt'
OPEN(11,FILE=fichero,status='old')
read(11,*)
read(11,*)datos
!find out number of data
ndat1=0
open(13,file=datos,status='old')
do 
  read(13,*,iostat=io)
  if(io/=0) exit
  ndat1=ndat1+1
enddo  
print *,'number of records',ndat1
close(13)

read(11,*)
read(11,*)genea
read(11,*)
read(11,*)tipomodelo

! Numero de efectos, numero de grupos geneticos (al menos uno), numero de caracteres, 
! Se asume que el ultimo caracter es el umbral
!ev-mad: number of efects (nfac) =number of covariates (ncov) + number of crossclassified effects (ncrossclassified)


READ(11,*)NFAC
READ(11,*)NYEAR
READ(11,*)NSEGM

!ev-mad ncov=number of covariates in the model
READ(11,*)NCOV
print *,'number of covariables: ',ncov
!
READ(11,*)NGRUP
if ((tipomodelo=='sire').and.(ngrup>0)) then
  print *,'Sire model and genetic groups are not allowed'
  stop
endif  
if ((tipomodelo=='g_usr').and.(ngrup>0)) then
  print *,'G_usr file model and genetic groups are not allowed'
  stop
endif  
READ(11,*)NCAR
print *,'number of observed traits',ncar
READ(11,*)ncarumb

IF(NCARumb.GT.nCAR)THEN
  PRINT *,'more threshold traits than total traits',NCARumb,ncar
  STOP
END IF
ncont=ncar-ncarumb
nres=0

allocate( nthr(ncarumb) )
nthr=0

!  numero de categorias (orden:primero 3 o mas, segundo 2)
maxthr=0
READ(11,*)(nthr(i),i=1,ncarumb)
do i=1,ncarumb
  if (nthr(i).eq.2) nres=nres+1
enddo 
  maxthr=maxval( nthr )

allocate( thr(ncar,maxthr),mint(ncar,maxthr) )
allocate( maxt(ncar,maxthr))

print *,'number of traits with var(e) constrained to 1 -->',nres

! numero de efectos permanentes
read(11,*) nper
! numero de efectos animales (correlacionados)
read(11,*) nefani
allocate( ifac(nfac+1)); ifac=0

! Numero de niveles de cada efecto (en los animales NO incluye grupo genetico)
READ(11,*)(IFAC(I),I=1,NFAC)
! Animales son los \FAltimos efectos 
! Permanente ("fijo aleatorio") los penultimos
NANIM=IFAC(NFAC)
do i=1,nefani
  ifac(nfac-nefani+i)=ifac(nfac-nefani+i)+ngrup
enddo

! Prepare address array (ifac)
NROW=0

! Ahora ifac(i) tiene la posicion del efecto (i-1)-esimo en su ultimo nivel
DO  I=1,NFAC
!ev-mad Modified to account for the levels of each covariate (1)
  if (i .le. ncov) then
    MM=1+nrow
    ifac(i)=nrow
    nrow=mm
  else
!
    MM=NROW+IFAC(I)
    IFAC(I)=NROW
    NROW=MM
  endif
enddo

! Numero total de ecuaciones
! USAR NROW PARA ALLOCATAR
IFAC(NFAC+1)=NROW
NANGRU=NANIM+NGRUP
PRINT *,'Animals, unknown parent groups =',NANIM,NGRUP
MAMI=0

do i=1,ncar
do iefani=1,nefani
write(aux,'(i2)')i
write(aux2,'(i2)')iefani
fichero='trends_trait_'//trim(adjustl(aux))//'_'//adjustl(aux2)
open(100+i*10+iefani,file=fichero)
enddo
enddo

allocate( b(nrow,ncar), &
     nzz(nrow),nfix(nfac,ndat1),y(ndat1,nCAR),MFAC(nfac), &
     vara(ncar*nefani,ncar*nefani), &
     varp(ncar,ncar,nper), &
     vare(nCAR,nCAR),vark(ncar) &
     ,se(nCAR,nCAR),ve(nCAR,nCAR) &
     ,sa(ncar*nefani,ncar*nefani),va(ncar*nefani,ncar*nefani)&
     ,reg(ndat1,ncar),IND(nCAR,nfac),KK2(nrow,nCAR) )

! Outputs
allocate( sumvara(ncar*nefani,ncar*nefani), &
     ssvara(ncar*nefani,ncar*nefani) )
allocate( sumvarp(ncar,ncar,nper), &
       ssvarp(ncar,ncar,nper) )
allocate( sumvare(ncar,ncar),ssvare(ncar,ncar) )

allocate( cora(ncar*nefani,ncar*nefani), &
     corp(ncar,ncar,nper) )
allocate( core(ncar,ncar) )
allocate( sumcora(ncar*nefani,ncar*nefani),sscora(ncar*nefani,ncar*nefani) )
allocate( sumcorp(ncar,ncar,nper), &
        sscorp(ncar,ncar,nper) )
allocate( sumcore(ncar,ncar),sscore(ncar,ncar) )
allocate( sumb(nrow,ncar),ssb(nrow,ncar) )
allocate( h2(ncar*nefani)      )
allocate( posvara1(nefani,ncar),posvara2(nefani,ncar) )
allocate( iped(nanim,5))
allocate( xmend(nanim,ncar,nefani))

allocate(sumy(ncar,nefani,nsegm+1,nyear))
allocate(bv(nanim,ncar,nefani,nsegm+1))
allocate(ssumyeee(ncar,nefani,nsegm+1,nyear))
allocate(sssumyeee(ncar,nefani,nsegm+1,nyear))
allocate(ny(nyear))
allocate(varyear(ncar,nefani,nyear))
allocate(svaryear(ncar,nefani,nyear))
allocate(ssvaryear(ncar,nefani,nyear))
allocate(ssvar(ncar,nefani,nyear))
allocate(ndesc(nanim))
allocate(nwei(nyear))
allocate(nsel(nyear))
allocate(sumwei(ncar,nefani,nyear))
allocate(sumsel(ncar,nefani,nyear))
allocate(differ(ncar,nefani,nyear))
allocate(intens(ncar,nefani,nyear))
allocate(sdiffer(ncar,nefani,nyear))
allocate(sintens(ncar,nefani,nyear))
allocate(ssdiffer(ncar,nefani,nyear))
allocate(ssintens(ncar,nefani,nyear))


! allocate linked list
allocate(zhz(maxrec),ifirst(nrow),ivcol(0:maxrec), &
     inext(0:maxrec) )

allocate( dia(nrow),yy(nrow,ncar) )
!EV
allocate(d(ncar), v(ncar,ncar))
!ev_mad para la cov
allocate(valor(ndat1,nfac))
valor=0d0

ndesc=0
ssvar=0
ssvaryear=0
svaryear=0
varyear=0

!     Set to zero (f90 style)
sumvara=0.d0
ssvara=0.d0
sumvarp=0.d0
ssvarp=0.d0
sumvare=0.d0
ssvare=0.d0
sumcora=0.d0
sscora=0.d0
sumcorp=0.d0
sscorp=0.d0
sumcore=0.d0
sscore=0.d0
sumb=0.d0
ssb=0.d0
d=0.d0;v=0.d0
var=0.d0
se=0.d0
ve=0.d0
sa=0.d0
va=0.d0
cora=0.d0
corp=0.d0
core=0.d0
h2=0.d0
b=0.d0
posvara1=0
posvara2=0
DIA=0.E+00
YY=0.E+00
reg=0.d0
y=0.d0
nfix=0
nzz=0.d0
IFIRST=0
ssumyeee=0
sssumyeee=0


!EV




! different models for each trait
DO I=1,NCAR
  READ(11,*)(IND(I,J),J=1,NFAC)
  print *, 'Model for trait',i
  PRINT *,(IND(I,J),J=1,NFAC)
ENDDO
! Reading task
read(11,*)
read(11,*) task !VCE or BLUP
if(task=='BLUP')VCE=.false.
if (VCE) print *,'Estimating variance components'
if (.not.(VCE)) print *,'Variance components fixed'

! Reading features of the Gibbs sampling
read(11,*)
read(11,*) imue
read(11,*)
read(11,*) lap
read(11,*)
read(11,*) icad
imue=imue/icad
lap=lap/icad
print *,'Total n of iterations:',imue*icad
print *,'Burn-in: discarded for results.txt and solutions.txt',lap*icad
print *,'Thin interval',icad
print *,'Total n of samples in the Gibbs sampler:',imue

! Reading initial variances
! Genetic
read(11,*)
do i=1,ncar*nefani
  read(11,*) vara(i,1:ncar*nefani) 
enddo
! call pos_def(vara(1:ncar*nefani,1:ncar*nefani),'vara not semi-positive def',1.d-20,stat)        
! if (stat) stop
! Random environmental (permanent)
read(11,*)
do iper=1,nper
  read(11,*)
  do i=1,ncar
    read(11,*) varp(i,1:ncar,iper)
  enddo
!  call pos_def(varp(1:ncar,1:ncar,iper),'varp not semi-positive def',1.d-20,stat)        
!  if (stat) stop
enddo
! Residual
read(11,*)
do i=1,ncar
  read(11,*) vare(i,1:ncar)
enddo


! call pos_def(vare(1:ncar,1:ncar),'vare not semi-positive def',1.d-20,stat)        
! if (stat) stop

! reading initial, optional random seeds
call getoption_unit(11,'RandomSeeds',pos1,x=seeds)
if(pos1/=-1) then
	if(pos1==3)then
		print *,'putting random seeds',int(seeds)
		call init_seeds(int(seeds(1)),int(seeds(2)),int(seeds(3)))
	endif
endif


close(11)

!-------------------------------------
! End reading parameter file
!-------------------------------------

!EV Continuaci\F3n de un proceso inacabado
!AL removed 5/4/06
!print*,&
!'\BFInicio del proceso (1),continuacion de un proceso inacabado(2) o continuacion de un proceso finalizado (3)?'
!read (*,*) respuesta
respuesta=1

!---------------
!-> Preliminares
!---------------



! funcion para la posicion de las varianzas
do iefani=1,nefani
  do j=1,ncar
    posvara1(iefani,j)=(iefani-1)*ncar+j
    posvara2(iefani,j)=(iefani-1)*ncar+j
  enddo
enddo

!EVcall crearotulo(ncar,nefani,nper,rotulo,pos1)
!EVwrite(15,'(100a15)')(rotulo(i),i=1,pos1)

!--------------------------------
! --> Loop Lectura genealogia <--
!--------------------------------

! El fichero debe estar ordenado por animal para chequeos

!     SE LEE FICHERO GEN Y SE CONSTRUYE LA INVERSA DE A
!       SIN CONSANGUINIDAD Y GPD

! --> Modelo animal (incluye grupos geneticos)
if (tipomodelo.eq.'animal') then 
  ny=0
  print *,'animal model with genetic groups '
  OPEN(12,FILE=GENEA,status='old')
  3 READ(12,*,END=4)IA,IP,IM,ig,iy
    iped(ia,1)=ia
    iped(ia,2)=ip
    iped(ia,3)=im
    iped(ia,4)=ig
    iped(ia,5)=iy
    ny(iy)=ny(iy)+1
    if (ip.le.nanim) then
       ndesc(ip)=ndesc(ip)+1
    endif
    if (im.le.nanim) then
       ndesc(im)=ndesc(im)+1
    endif

! VERIFICACION DE LA ESTRUCTURA DE LOS DATOS
  MAMI=MAMI+1
  if(mod(mami,10000)==0) print *,'pedigree file, record',mami


  IF((IA.GT.NANIM).OR.(IP.GT.NANGRU).OR.(IM.GT.NANGRU).OR. &
    (IA.NE.MAMI).OR.(IP.LT.1).OR.(IM.LT.1))THEN
    PRINT *,'non-sorted or incorrect genealogy file?',IA,IP,IM
    STOP
  ENDIF

! Almacenamos A-1 repetidamente
  do i=1,nefani
    ! Si hay genealogia desconocida 
    XM=0.d0
    IF(IP.GT.NANIM) XM=XM+1.
    IF(IM.GT.NANIM) XM=XM+1.
    X=4./(XM+2.)
    ! Posiciones
    IA1=IA+IFAC(NFAC-nefani+i)
    IP1=IP+IFAC(NFAC-nefani+i)
    IM1=IM+IFAC(NFAC-nefani+i)

    ! Almacenamos el valor del animal en la diagonal (es la diagonal de A-1)
    DIA(IA1)=DIA(IA1)+X

    ! Animal - Padre y madre
    X2=X/2.
      CALL LINKS(IP1,IA1,-X2)
      CALL LINKS(IA1,IP1,-X2)
      CALL LINKS(IM1,IA1,-X2)
      CALL LINKS(IA1,IM1,-X2)

    ! Padre-madre
    X2=X/4.
    DIA(IP1)=DIA(IP1)+X2
    DIA(IM1)=DIA(IM1)+X2
    IF(IM1/=IP1)THEN
      CALL LINKS(IM1,IP1,X2)
      CALL LINKS(IP1,IM1,X2)
    ELSE
      DIA(IM1)=DIA(IM1)+2*X2
    ENDIF
  enddo
  GOTO 3
4 CLOSE(12)
endif
nwei=0
nsel=0
do i=1,nanim
 nwei(iped(i,5))=nwei(iped(i,5))+ndesc(i)
 if (ndesc(i).ge.1) then
 nsel(iped(i,5))=nsel(iped(i,5))+1
 endif
enddo    

    
    
! --> sire model without genetic groups 

if(tipomodelo.eq.'sire') then 
  print*, 'sire model'
  OPEN(12,FILE=GENEA,status='old')
55  READ(12,*,END=66)IS,ISS,IMGS
  ! VERIFICACION DE LA ESTRUCTURA DE LOS DATOS
  mami=MAMI+1
  if(mod(mami,10000)==0) print *,'pedigree file, record',mami
  IF((IS.GT.NANIM).OR.(ISS.GT.NANGRU).OR.(IMGS.GT.NANGRU).OR. &
   (IS.NE.MAMI) .OR.(ISS.GT.NANIM) .OR. (IMGS.GT.NANIM)) THEN
    PRINT *,'non-sorted or incorrect genealogy file?A',IS,ISS,IMGS
    STOP
  ENDIF
    
  ! Almacenamos A-1 repetidamente
  do i=1,nefani
    ! Posiciones
    IS1=IS+IFAC(NFAC-nefani+i)
    ISS1=ISS+IFAC(NFAC-nefani+i)
    IMGS1=IMGS+IFAC(NFAC-nefani+i)

    ! si sire y MGS known
    if ((ISS .ne. 0) .and. (IMGS .ne. 0)) THEN
      X=16./11.
      ! Almacenamos el valor del SIRE en la diagonal(es la diagonal de A-1)
      DIA(IS1)=DIA(IS1)+X
      ! (SS,S) o (S,SS)
      CALL LINKS(ISS1,IS1,-X/2.)
      CALL LINKS(IS1,ISS1,-X/2.)
      ! (MGS,S) o (S,MGS)
        CALL LINKS(IMGS1,IS1,-X/4.)
        CALL LINKS(IS1,IMGS1,-X/4.)
       ! (SS,SS)
      DIA(ISS1)=DIA(ISS1)+X/4.
      ! (SS,MGS) O (MGS, SS)
      IF(ISS1/=IMGS1)THEN
        CALL LINKS(ISS1,IMGS1,X/8.)
        CALL LINKS(IMGS1,ISS1,X/8.)
      END IF
      IF (ISS1==IMGS1)THEN
        DIA(ISS1)=DIA(ISS1)+X/8.
        DIA(IMGS1)=DIA(IMGS1)+X/8.
      END IF
      ! (MGS,MGS)
      DIA(IMGS1)=DIA(IMGS1)+X/16.
    END IF

    ! si solo sire OF SIRE conocido:
    if ((ISS .ne. 0 ) .and. (IMGS .eq. 0)) THEN
      ! Almacenamos el valor del SIRE en la diagonal
      ! (es la diagonal de A-1)
      X=4./3.
      DIA(IS1)=DIA(IS1)+X
      ! (SS,S) o (S,SS)
      CALL LINKS(ISS1,IS1,-X/2.)
      CALL LINKS(IS1,ISS1,-X/2.)
      ! (SS,SS)
      DIA(ISS1)=DIA(ISS1)+X/4.
    END IF
    ! si solo mgs conocido:
    if ((ISS .eq. 0 ) .and. (IMGS .ne. 0)) THEN
      X=16./15.
      DIA(IMGS1)=DIA(IMGS1)+X
      ! (MGS,S) o (S,MGS)
        CALL LINKS(IMGS1,IS1,-X/4.)
        CALL LINKS(IS1,IMGS1,-X/4.)
      ! (MGS,MGS)
      DIA(IMGS1)=DIA(IMGS1)+X/16.
    END IF 
    ! si ninguno conocido
    if ((ISS .eq. 0) .and. (IMGS .eq. 0)) THEN
      X=1.
      DIA(IS1)=DIA(IS1)+X
    END IF
    ENDDO
  GOTO 55
66 CLOSE(12)
endif

! --> user_defined
if(tipomodelo=='g_usr') then
    mami=0
    print *,'user defined G inverse'
    ! this file has to be lower or upper stored or mixed
    ! but not both !!
    OPEN(12,file=genea,status='old')
155  READ(12,*,END=166)IA,IP,x ! i,j,a values
  mami=mami+1
  if(mod(mami,10000)==0) print *,'g_usr file, record',mami
    
  ! Almacenamos A-1 repetidamente
  do i=1,nefani
    ! Posiciones
    IA1=IA+IFAC(NFAC-nefani+i)
    IP1=IP+IFAC(NFAC-nefani+i)

    if(ia1/=ip1) then
    	call links(ip1,ia1,x)
    	call links(ia1,ip1,x)
    else
    	DIA(IA1)=DIA(IA1)+X
    endif

    ENDDO
  GOTO 155
166 CLOSE(12)
endif

IF( (tipomodelo/='g_usr') .and. (NANIM.NE.MAMI) )THEN
  PRINT *,'different number of animals in genealogy than in par file',MAMI,NANIM
  STOP
ENDIF


!----------------------------------
! -- Fin Loop Lectura genealogia --
!----------------------------------


!-------------------------------------
! --> Lectura del fichero de datos <--
!-------------------------------------

!     LECTURA DEL FICHERO DE DATOS
OPEN(13,FILE=datos,status='old')
NDAT1=0
5 NDAT1=NDAT1+1
  if(mod(ndat1,10000)==0) print *,'Data file, record',ndat1


  !  Lectura de datos: covariables(mod_ev_mad), efectos fijos, fijo aleatorio,animal, caracteres
  !  leo todo en mfac (a diferencia de antes q eran dos vectores, ia y mfac)
  !      READ(13,*,END=6)(MFAC(I),I=nfac-nefani+1,NFAC),
  !       animales
  !     +(MFAC(I),I=1,NFAC-nefani)
  !       fijos y permanentes
  !     1,(REG(NDAT1,I),i=1,ncar)
!ev_mad the vector valor takes the value of the covariates in their positions
!mfac vector takes the value of each crossclassified effect
  READ(13,*,END=6) (valor(ndat1,i),i=1,ncov),(MFAC(I),I=ncov+1,NFAC), & !efectos
            (REG(NDAT1,I),i=1,ncar)  !registros


  do i=nfac-nefani+1,nfac
    IF((mfac(i).LT.1).OR.(mfac(i).GT.NANIM))THEN
     PRINT *,'animal=0 or higher number than expected',mfac(i),NANIM
     STOP                       
    END IF
  enddo
!ev_mad: I modified it in order to skip the covariates
  DO  I=1+NCOV,NFAC-nefani
    IF((MFAC(I).LT.1).OR.(MFAC(I).GT.IFAC(I+1)-IFAC(I)))THEN
      PRINT *,'level of effect',I,' beyond bounds in record',  &
       ndat1,' : ',IFAC(I+1)-IFAC(I),' --->',MFAC(I)
      STOP
    END IF
  enddo
!ev_mad
!valor vector takes the value of the covariate in their corresponding positions and 1 in the corresponding positions of
!the crossclassified effects
!mfac vector takes the value of 1 in the positions corresponding to the covariates and the level of each crossclassified
!effect in ther corresponding positions for the crossclassified effects
!Example: data 1: 278 250 5 (cov1 cov2 ef1), valor will contain (270 250 1), and mfac (1 1 5)


        do j=1,nfac
                if (j .le. ncov) then
                mfac(j)=1
                else
                valor(ndat1,j)=1.d0
                endif
        enddo


  ! Nfix almacena toda la estructura de datos para todos los registros. O sea [X;Z]
  ! En realidad es [Z;X] modificarlo a [Z1:Zn;X]
  ! loop a lo largo de todos los ef animales
  do i=1,nefani
    NFIX(i,NDAT1)=mfac(nfac-nefani+i)
  enddo
  DO  I=1,NFAC-nefani
    NFIX(I+nefani,NDAT1)=MFAC(I)
  enddo
  
  ! SE CONSTRUYE RHS: XX, ZX 

  do i=1,nfac
    pos1=ifac(i)+mfac(i)
    ! ef fijo/permanente -> X'X
!ev_mad *valor(ndat1,i) was added. for covariates, we multiply by the value of the covariate**2,that is, valor*valor)
    if (i.le.(nfac-nefani)) dia(pos1)=dia(pos1)+valor(ndat1,i)*valor(ndat1,i)
        ! ef aditivo -> Z'Z
    if (i.gt.(nfac-nefani)) nzz(pos1)=nzz(pos1)+valor(ndat1,i)*valor(ndat1,i)
    do j=1,nfac
      pos2=ifac(j)+mfac(j)
      ! -> X'Z, Z1'Z2,etc
          !ev_mad valor(ndat1,i)*valor(ndat1,j) was added. If pos1 corresponds to a covariate and pos2 corresponds to another effect 
          !(not covariate), the value of the covariate*1 is added. If pos1 and pos2 are not covariates, 1*1 is added
          if (pos1/=pos2) CALL LINKS(pos1,pos2,(valor(ndat1,i)*valor(ndat1,j)))
        enddo
  enddo

GOTO 5
6 CLOSE(13)
NDAT1=NDAT1-1

!---------------------------------------
! -- Fin Lectura del fichero de datos --
!---------------------------------------


       PRINT*,'number of records, non-null elements = ',NDAT1,NPLACE*2

!===========================================================================
!
!   B.-GIBBS


! \C1 couple of preliminaries...

! kk2 es un indicador de si el muestreo en 
!  la posicion (efecto, caracter) es v\E1lida (.true.) o no (.false.)

      KK2=.false.

      DO I=1,NFAC
        DO J=1,NCAR
          IF (IND(J,I).EQ.1) THEN
            DO K=IFAC(I)+1,IFAC(I+1)
               KK2(k,j)=.true.
            ENDDO
          ENDIF
        ENDDO
      ENDDO

       nanim=nanim+ngrup
!      Posicion de partida de los animales en las MME
      ilev=ifac(nfac-nefani+1)


! Varianzas iniciales


!- Inicializo los umbrales
! Primer umbral es 0, segundo (si >=2 categorias) 1, y ultimo cuasi-infinito 
! loop para todos los car. umbrales
      do j=1,ncarumb
        do i=1,nthr(j)-1
          thr(j,i)=1.d0*(i-1)
        enddo
        thr(j,nthr(j))=999.d0
	!AL 6/04/06 to avoid jumpins in binary traits
	if (nthr(j)==2) thr(j,nthr(j))=liabilitybound !for binary traits only!!!

        mint(j,1)=-999.d0
        maxt(j,1)=0.d0
        if(nthr(j).gt.2) then
          mint(j,2)=0.d0
!          maxt(j,2)=10.d0
          maxt(j,2)=1.d0
        endif
        maxt(j,nthr(j))=thr(j,nthr(j))

!  si >2 categorias escapo los 2 primeros y el ultimo;
! si no, el 1o y el ultimo
        if (nthr(j).gt.2) then
          do i=3,nthr(j)-1
            mint(j,i)=thr(j,i-1)+1.d-2
            maxt(j,i)=thr(j,i)-1.d-2
          enddo
        else
          do i=2,nthr(j)-1
            mint(j,i)=thr(j,i-1)+1.d-2
            maxt(j,i)=thr(j,i)-1.d-2
          enddo
        endif
        mint(j,nthr(j))=thr(j,nthr(j)-1)+1.d-2
      enddo


!---------------------------------------------------------------------
!EV Continuaci\F3n de un proceso inacabado
!print*,&
!'\BFInicio del proceso (1),continuacion de un proceso inacabado(2) o continuacion de un proceso finalizado (3)?'
!read (*,*) respuesta
if (respuesta.eq.1) then
  close(15)
  close(16)
  close(20)
  open(15,file='samples.txt',status='replace')
  open(16,file='thresholds.txt',status='replace')
  open(20,file='samplesFE.txt',status='replace')
  call crearotulo(ncar,nefani,nper,rotulo,pos1,tipomodelo)
  write(15,'(100a15)')(rotulo(i),i=1,pos1)
endif
if (respuesta.eq.3) then
  print*,' \BFcuantas iteraciones mas?'
  read*, itermas
endif
if ((respuesta.eq.2).or.(respuesta.eq.3)) then 
  open(21, file ='continuacion.txt')
  read(21,'(21x,i10)')itercor
  read(21,*)
  read(21,*)
  read(21,*)
  read(21,*)  
  read(21,*)
    do i=1, ncar*nefani
      read(21,'(20f15.8)')(vara(i,j),j=1,ncar*nefani)
    enddo
    do k=1,nper
      read(21,*)
	    do i=1,ncar
		  read(21,'(20f15.8)')(varp(i,j,k),j=1,ncar)
	    enddo
    enddo
      read(21,*)
        do i=1,ncar
		  read(21,'(20f15.8)') (vare(i,j),j=1,ncar)
        enddo
      read(21,*)
        do j=1,ncarumb
		  read(21,'(20f15.8)') (thr(j,i),i=1,nthr(j))
		  print*, (thr(j,i),i=1,nthr(j))
        enddo
       read(21,'(8x,i10)')  x1
       print*, x1
	    read(21,*)
          do i=1, nrow
		   read(21,*)b(i,1:ncar)
          enddo
endif
close(21)



!----------------------------------------------------------------------

!     Muestreo de Gibbs
if (respuesta.eq.1) itercor=0
if (respuesta.eq.3) imue=itermas+itercor

!EV
if((respuesta.eq.2).or.(respuesta==3)) then 
  rewind(15)
  read(15,'(a10)')kkk
!  print*, 'kk',kkk
    do i=1, itercor/icad
      read(15,'(a10)')kkk
!      print*,kkk
!     pause
    enddo
endif


      do 9999 ijk=(itercor/icad)+1,imue
      do 9998 kji=1,icad
!      pause
	print *,ijk,kji

! ----- Invierto matrices de covarianzas -----------------

        call ginv1(vare,ncar,ncar,tol,irango)
        call ginv1(vara,ncar*nefani,ncar*nefani,tol,irango)

        do i=1,nper
          se=varp(:,:,i)
          call ginv1(se,ncar,ncar,tol,irango)
          varp(:,:,i)=se
        enddo

! ------- aumento de datos para datos faltantes o censurados ---------------
!         y para caracteres discontinuos


!  Caracteres lineales
! Es un poco embarullado y me pierdo un poco
        do j=1,ncont
        do i=1,ndat1

!EV inicializo vark
!         vark=0.d0
         if (reg(i,j).lt.0.01) then
              xm=0.
              do k=1,ncar
! vark es una variable de trabajo
               vark(k)=y(i,k)
              enddo

!              do jk=1,nfac-1
              do jk=1,nfac-nefani
! b es el vector de soluciones
! por tanto esto es Xb
!ev_mad he anhadido *valor(i,jk)
                  xm=xm+b(nfix(jk+nefani,i)+ifac(jk),j)*valor(i,jk)
                  do k=1,ncar
                    if (k.ne.j) then
                      vark(k)=vark(k)-b(nfix(jk+nefani,i)+ifac(jk),k)*valor(i,jk)
                    endif
                  enddo
              enddo
! y esto es Zb (son los animales)
! loop para efectos maternos 
!ev_mad he anhadido *valor(i,nfac-nefani+jk)
              do jk=1,nefani
                xm=xm+b(nfix(jk,i)+ifac(nfac-nefani+jk),j)*valor(i,nfac-nefani+jk)
!              xm=xm+b(nfix(1,i)+ifac(nfac),j)
                do k=1,ncar
                  if (k.ne.j) then
!ev_mad he anhadido *valor(i,nfac-nefani+jk)
                    vark(k)=vark(k)-b(nfix(jk,i)+ifac(nfac-nefani+jk),k)*valor(i,nfac-nefani+jk)
!                    vark(k)=vark(k)-b(nfix(1,i)+ifac(nfac),k)
                  endif
                enddo
            enddo
! e~MVN(e,R)
              do k=1,ncar
                if (k.ne.j) then
                  xm=xm-vark(k)*vare(k,j)/vare(j,j)
                endif
              enddo
              var=1.d0/vare(j,j)

!este para missing value
			  if (reg(i,j).gt. -0.01) then
				  call normal(x1,u)
				  y(i,j)=xm+u*sqrt(var)
			  elseif (reg(i,j) .lt. -0.01) then
!EV  para datos censurados, se codifica la cota inferior como -cota
				  cota1=-reg(i,j) -xm
				  cota2=999.d0 -xm
				  cota1=cota1/sqrt(var)
				  cota2=cota2/sqrt(var)
				  call gen_trunc_normal(cota1,cota2,u,x1)
				  y(i,j)=xm+u*sqrt(var)
			  endif
! Var usa una formulacion equivalente q se basa en las 
! propiedades de la inversa de MVN.
         else
              y(i,j)=reg(i,j)
         endif
        enddo
        enddo

! -->   aumento de los caracteres umbrales (ultimos)

        do j=1,ncarumb 
          poscar=ncont+j
          do i=1,ndat1

              xm=0.
              do k=1,ncar
               vark(k)=y(i,k)
              enddo

! Xb
!              do jk=1,nfac-1
              do jk=1,nfac-nefani
!ev_mad He anhadido *valor(i,jk)
                  xm=xm+b(nfix(jk+nefani,i)+ifac(jk),poscar)*valor(i,jk)
                                  
                  do k=1,ncar
                    if (k.ne.poscar) then
!ev_mad He anhadido *valor(i,jk)
                      vark(k)=vark(k)-b(nfix(jk+nefani,i)+ifac(jk),k)*valor(i,jk)
                    endif
                  enddo
              enddo
! Zu
              do jk=1,nefani
!ev_mad He anhadido *valor(i,nfac-nefani+jk)
                xm=xm+b(nfix(jk,i)+ifac(nfac-nefani+jk),poscar)*valor(i,nfac-nefani+jk)
                do k=1,ncar
                  if (k.ne.poscar) then
!ev_mad He anhadido *valor(i,nfac-nefani+jk)
                    vark(k)=vark(k)-b(nfix(jk,i)+ifac(nfac-nefani+jk),k)*valor(i,nfac-nefani+jk)
                                  endif
                enddo
            enddo
! e~MVN(e,R)
              do k=1,ncar
                if (k.ne.poscar) then
                  xm=xm-vark(k)*vare(k,poscar)/vare(poscar,poscar)
                endif
              enddo

              var=1.d0/vare(poscar,poscar)
          
!       Bajo que umbral esta?

              presenthr=reg(i,poscar)
! --> generacion de liability si 'missing value' : no hay limites, es igual que un
!  dato lineal
! ahora lo limito entre -999 y +999 (umbrales max y min)
              if (dble(presenthr).lt.0.01d0) then !missing or censored
	        if (dble(presenthr).gt.-0.01d0) then
		  !real missing value, no lower bound
                  cota1=(-999.d0-xm)
                  cota2=(999.d0-xm)
		  !AL 6/4/06 For binary traits to avoid going to infinite the 
  		  !upper and lower bound are 4 (4 s.d. from the threshold)
	          if (nthr(j)==2) then
		    cota1=-liabilitybound-xm
		    cota2=liabilitybound-xm	
		  endif
                  cota1=cota1/sqrt(var)
                  cota2=cota2/sqrt(var)
                  call gen_trunc_normal(cota1,cota2,u,x1)
!                 call normal(x1,u)
                  y(i,poscar)=xm+u*sqrt(var)
	        else !censored categorical trait - has to be more than 3 categories, otherwise censoring is meaningless
                  if (presenthr.eq.-1) then
                    cota1=-999.d0-xm
                    !AL 6/4/06 For binary traits to avoid going to infinite
		    if (nthr(j)==2) then
		    cota1=-liabilitybound-xm
		  endif
                  cota2=999.d0-xm
                  else
                    cota1=(thr(j,abs(presenthr)-1)-xm)
                    cota2=(999.d0-xm)
		    !AL 6/4/06 For binary traits to avoid going to infinite
		    if (nthr(j)==2) then
		      cota2=+liabilitybound-xm
		    endif
                  endif
                  cota1=cota1/sqrt(var)
                  cota2=cota2/sqrt(var)
                  call gen_trunc_normal(cota1,cota2,u,x1)
                  y(i,poscar)=xm+u*sqrt(var)
	        endif
		
		  
! --> generacion de liability cuando si se tiene el fenotipo
              else 
                if (presenthr.eq.1) then
                  cota1=-999.d0-xm
		  !AL 6/4/06 For binary traits to avoid going to infinite
		  if (nthr(j)==2) then
		    cota1=-liabilitybound-xm
		  endif
                  cota2=thr(j,presenthr)-xm
                else
                  cota1=(thr(j,presenthr-1)-xm)
                  cota2=(thr(j,presenthr)-xm)
		  !AL 6/4/06 For binary traits to avoid going to infinite
		  if (nthr(j)==2) then
		    cota2=+liabilitybound-xm
		  endif
                endif
                cota1=cota1/sqrt(var)
                cota2=cota2/sqrt(var)
                call gen_trunc_normal(cota1,cota2,u,x1)
                y(i,poscar)=xm+u*sqrt(var)

!-- Almaceno minimos y maximos para el muestreo de los umbrales
!    solo cuando son datos "verdaderos"
                if (y(i,poscar).lt.mint(j,presenthr)) then
                  mint(j,presenthr)=y(i,poscar)
                endif
                if (y(i,poscar).gt.maxt(j,presenthr)) then
                  maxt(j,presenthr)=y(i,poscar)
                endif
                if((y(i,poscar).lt.0).and.(presenthr.eq.2)) then
                  print *,y(i,poscar),reg(i,poscar),cota1,cota2,u,xm,var
                  print *,'kk1'
                  print *,vare(poscar,poscar)
                  stop
                endif
              endif
          enddo
        enddo
!        Fin del loop de aumento de datos

!
! --------- Muestreo de umbrales ------------------

! Escapo el primer umbral porque es 0 y el ultimo pq es +infinito
! Escapo los dos primeros umbrales (0,1) si hay > 2 categorias

        do j=1,ncarumb 
! Para el umbral 1 (=0)
          maxt(j,1)=0.d0              
! Para el umbral 2 (=10)
          if(nthr(j).gt.2) then
            mint(j,2)=0.d0
            maxt(j,2)=1.d0      
          endif
          if (nthr(j).gt.2) then
            do i=3,nthr(j)-1
              cota1=maxt(j,i)
              cota2=mint(j,i+1)
              if (cota2.lt.cota1) then
                print *,cota1,cota2,'cotas'
                stop
              endif
              call unif(x1,u)
              thr(j,i)=cota1+u*(cota2-cota1)
            enddo

! Borro minimos y maximos para el muestreo de los umbrales
            do i=1,nthr(j)
              maxt(j,i)=-1000.d0
              mint(j,i)=+1000.d0
            enddo
            mint(j,1)=-999.d0
            maxt(j,nthr(j))=999.d0
          endif
        enddo

!	print *,'pasa2',ssumyeee(1,1,1,1)

! ------------- contruccion de rhs --------------------------------

        yy=0.d0

        do i=1,ndat1
           do k=1,ncar
               do j=1,nfac-nefani
! Construye X'y en una columna por caracter
!ev_mad He anhadido valor(i,j)
                 yy(nfix(j+nefani,i)+ifac(J),k)=  &
                 yy(nfix(j+nefani,i)+ifac(J),k)+y(i,k)*valor(i,j)
               enddo               
! Contruye Z'y id.
!               yy(nfix(1,i)+ilev,k)=yy(nfix(1,i)+ilev,k)+y(i,k)
               do j=1,nefani
                 pos1=ifac(nfac-nefani+j)
!ev_mad He anhadido valor(i,j+nfac-nefani)
                 yy(nfix(j,i)+pos1,k)=yy(nfix(j,i)+pos1,k)+y(i,k)*valor(i,j+nfac-nefani)
               enddo
           enddo
        enddo

!	print *,'pasa3',ssumyeee(1,1,1,1)
! ---------- muestreo de efectos fijos -----------------------------------------

         do 1005 i=1,ifac(nfac-nper-nefani+1)
           do 5001 j=1,ncar
               b(i,j)=0.d0
! Primero asigna a b el valor del RHS (xR-1y)
               do 2001 k=1,ncar
                 b(i,j)=b(i,j)+yy(i,k)*vare(j,k)
                 if (j.ne.k) then
! Y quita el efecto de correlacion de los otros caracteres
                   b(i,j)=b(i,j)-b(i,k)*vare(j,k)*dia(i)
                 endif
2001           continue

! Le resta las soluciones de los demas efectos q estan acumuladas en ZHZ
               iplace=ifirst(i)
1006           if(iplace.gt.0)then
                 ji=ivcol(iplace)               
                   do 2002 k=1,ncar
                     b(i,j)=b(i,j)-b(ji,k)*zhz(iplace)*vare(j,k)             
2002               continue
                 iplace=inext(iplace)
                 goto 1006
               endif
! Sistema para que muestree los efectos segun el modelo
               KK=KK2(I,J)
               IF (KK) THEN
                 b(i,j)=b(i,j)/(dia(i)*vare(j,j))
                 var=1./(dia(i)*vare(j,j))
                u=xnor(x1)
!               call normal(x1,u)
                b(i,j)=b(i,j)+u*sqrt(var)
               ELSE
                 B(I,J)=0.d0
               ENDIF
5001       continue
1005     continue
!	print *,'pasa4',ssumyeee(1,1,1,1)
!-------muestreo de efectos fijos aleatorios---------------------------
! loop a lo largo de todos los efectos permanentes
        do iper=1,nper
!        do i=ifac(nfac-nper+iper-1)+1,ifac(nfac-nper+iper)
!        do i=ifac(nfac-nper+iper-nefani)+1,ifac(nfac-nper+iper-nefani+1)
        pos1=ifac(nfac-nper+iper-nefani)+1
      pos2=ifac(nfac-nper+iper-nefani+1)
        do i=pos1,pos2
!           PRINT *,I,DIA(I)
          do j=1,ncar
             b(i,j)=0.d0
             do k=1,ncar
               b(i,j)=b(i,j)+yy(i,k)*vare(j,k)
               if (j.ne.k) then
                 b(i,j)=b(i,j)-b(i,k)*vare(j,k)*dia(i)
                 b(i,j)=b(i,j)-b(i,k)*varp(j,k,iper)
               endif
             enddo    
!            PRINT *,IFIRST(I)
             iplace=ifirst(i)
1426         if(iplace.gt.0)then
               ji=ivcol(iplace)
         
               do k=1,ncar
                 b(i,j)=b(i,j)-b(ji,k)*zhz(iplace)*vare(j,k)
               ENDDO
               iplace=inext(iplace)
               goto 1426
             endif
! Sistema para que muestree los efectos segun el modelo

             KK=KK2(I,J)
             IF(KK) then
               b(i,j)=b(i,j)/(dia(i)*vare(j,j)+varp(j,j,iper))
               var=1./(dia(i)*vare(j,j)+varp(j,j,iper))
               u=xnor(x1)
!              call normal(x1,u)
               b(i,j)=b(i,j)+u*sqrt(var)
             else
               b(i,j)=0.d0
             endif
        
          ENDDO
        ENDDO
        enddo
!	print *,'pasa5',ssumyeee(1,1,1,1)
!------Muestreo de aleatorios -------------------------------------------

       do iefani=1,nefani
       pos1=ifac(nfac-nefani+iefani)+1
       pos2=ifac(nfac-nefani+iefani+1)
       do 1007 i=pos1,pos2
         do 5002 j=1,ncar
!           posvara1=(iefani-1)*ncar+j
           b(i,j)=0.
           do 2003 k=1,ncar
             b(i,j)=b(i,j)+yy(i,k)*vare(j,k)
             if (j.ne.k) then
!              covarianzas residuales
               b(i,j)=b(i,j)-b(i,k)*vare(j,k)*nzz(i)
             endif
2003       continue
!         esto es la correlacion de los otros efectos aditivos
!         del mismo individuo
           do l=1,nefani
             do k=1,ncar
               if (posvara1(iefani,j).ne.posvara2(l,k)) then
!                 posicion de la solucion del animal
!                 b(i,j)=b(i,j)-b(nanim*(l-1)+i,k)*
!                 posani=i-ifac(nfac-nefani+iefani)+ifac(nfac-nefani+l)
                 b(i,j)=b(i,j)- &
                b(i-ifac(nfac-nefani+iefani)+ifac(nfac-nefani+l),k)* &
                vara(posvara1(iefani,j),posvara2(l,k))*dia(i) 
               endif
             enddo
           enddo

!     barrido del resto 
           iplace=ifirst(i)
1008       if(iplace.gt.0)then
             ji=ivcol(iplace)
             do k=1,ncar
! fijos / permanentes/Z1'Z2
!               if((ji.le.pos1).or.(ji.gt.pos2))then
               if((ji.lt.pos1).or.(ji.gt.pos2))then
                 b(i,j)=b(i,j)-b(ji,k)*zhz(iplace)*vare(j,k)
               else 
! A-1*G-1
                 do l=1,nefani
!     posicion de la solucion del animal
!                   b(i,j)=b(i,j)-b(nanim*(l-1)+ji,k)*zhz(iplace)*
!                 posani=ji-ifac(nfac-nefani+iefani)+ifac(nfac-nefani+l)
                   b(i,j)=b(i,j)- &
      b(ji-ifac(nfac-nefani+iefani)+ifac(nfac-nefani+l),k)*zhz(iplace)* &
                        vara(posvara1(iefani,j),posvara2(l,k))
                 enddo
               endif
             enddo
             iplace=inext(iplace)
             goto 1008
           endif
! muestreo segun modelo...
           KK=KK2(I,J)
           IF(KK) then
             var=1./(nzz(i)*vare(j,j)+dia(i)* &
             vara(posvara1(iefani,j),posvara1(iefani,j))) 
             b(i,j)=b(i,j)/(nzz(i)*vare(j,j)+ &
             dia(i)*vara(posvara1(iefani,j),posvara1(iefani,j)))
             u=xnor(x1)
!            call normal(x1,u)
             b(i,j)=b(i,j)+u*sqrt(var)
           else
             b(i,j)=0.d0
           endif

5002     continue
1007   continue
       

! Pone a 0 el ultimo grupo genetico
        do i=1,ncar 
          if(tipomodelo.eq.'animal')b(ifac(nfac-nefani+iefani+1),i)=0.d0
        enddo
       enddo
!        fin loop nefani

!-- SAMPLING VARIANCE COMPONENTS ONLY IF DESIRED --
       VarianceComponentsEstimated: if(VCE) then

!	print *,'pasa6',ssumyeee(1,1,1,1)
!------Muestreo de varianzas residuales----------------------------------


! Utiliza flat priors para todas las varianzas

        se=0.d0

! Calculo de residuos
        do 1011 i=1,ndat1
           do 2009 k=1,ncar
              vark(k)=y(i,k)
!  Fijos y permanentes
              do 1019 j=1,nfac-nefani
!ev_mad he anhadido valor(i, j)
                 vark(k)=vark(k)-b(nfix(j+nefani,i)+ifac(j),k)*valor(i,j)
1019          continue
!  Aleatorios
              do j=1,nefani
!ev_mad he anhadido valor(i, j+nfac-nefani)
                vark(k)=vark(k)-b(nfix(j,i)+ifac(nfac-nefani+j),k)*valor(i,j+nfac-nefani)
              enddo
!                vark(k)=vark(k)-b(nfix(1,i)+ifac(nfac),k)
2009      continue

!   Suma de cuadrados
          do 2010 k=1,ncar
            do 2010 l=1,ncar
             se(k,l)=se(k,l)+vark(k)*vark(l)
2010      continue
    
1011    continue


! Se ~ IW
		
        se=se/real(ndat1-NCAR-1)
!	  do
        if (nres.eq.0) then
        call ginv1(se,ncar,ncar,tol,irango)
        call wish(ncar,se,ve,ndat1-NCAR-1,x1)
        call ginv1(ve,ncar,ncar,tol,irango)
        else
!-- esta se usa para restringir la varianza residual del umbral a 1
          call inv_con_wish_multiple(ncar,nres,se,ve,ndat1-NCAR-1,x1)
        endif

!EV subrutina de denseop (Misztal)
!        call pos_def(ve(1:ncar,1:ncar),'ve no positiva def, corregida',1.d-6,stat)        
!        if (stat) npdve=npdve+1

!EV lo he metido dentro de un loop. Igual para el resto de las varianzas
!		call eigen(ve(1:ncar,1:ncar),d(1:ncar),v(1:ncar,1:ncar))
!		if (minval(d(1:ncar)).gt.1.d-1) exit
!		print*, 'eigenvalueve' , minval(d(1:ncar))
!		print*, 'matriz eig', ve(1:ncar,1:ncar), d(1:ncar)
!		pause
!	  enddo
!        unicaracter
!        call inv_con_wish(ncar,se,ve,ndat1-NCAR-1,x1)

        vare=ve
!	print *,'pasa7',ssumyeee(1,1,1,1)
	
!------Muestreo de varianzas permanentes----------------------------------
! loop para todos los permanentes
        do iper=1,nper
          se=0.d0

! Suma de cuadrados
!          do i=ifac(nfac-nper+iper-1)+1,ifac(nfac-nper+iper)
          pos1=ifac(nfac-nper+iper-nefani)+1
          pos2=ifac(nfac-nper+iper-nefani+1)
          do i=pos1,pos2
            do k=1,ncar
              do l=1,ncar
                 se(k,l)=se(k,l)+b(i,k)*b(i,l)
              enddo
            enddo
          enddo

! Muestreo
              se=se/real(pos2-pos1+1-NCAR-1)

!          do
            call ginv1(se,ncar,ncar,tol,irango)
            call wish(ncar,se,ve,pos2-pos1+1-NCAR-1,x1)
            call ginv1(ve,ncar,ncar,tol,irango)

!EV subrutina de denseop (Misztal)
!            call pos_def(ve(1:ncar,1:ncar),'vp no positiva def, corregida',1.d-6,stat)        
!            if (stat) npdvp=npdvp+1


!EV
!EV lo he metido dentro de un loop. Igual para el resto de las varianzas
!		    call eigen(ve(1:ncar,1:ncar),d(1:ncar),v(1:ncar,1:ncar))
!		    if (minval(d(1:ncar)).gt.1.d-1) exit
!			print*, 'eigenvaluevp' ,iper, minval(d(1:ncar))
!		    print*, 'matriz eig', ve(1:ncar,1:ncar), d(1:ncar)
!            pause
!		  enddo
          varp(:,:,iper)=ve
        enddo

!----------muestreo de varianzas aditivos-------------------------------
! 
        sa=0.d0

! Suma a'A-1a
! A-1 esta almacenado en ZHZ

!        trabajamos con A-1 abajo del todo
        do i=ifac(nfac)+1,ifac(nfac+1)
          do iefani=1,nefani
            pos1=i-ifac(nfac)+ifac(nfac-nefani+iefani)
! Elementos diagonales
            do iefani2=1,nefani
            pos2=i-ifac(nfac)+ifac(nfac-nefani+iefani2)
             do  k=1,ncar
!             posvara1=(iefani-1)*ncar+k
              do  l=1,ncar
!              posvara2=(iefani2-1)*ncar+l
                sa(posvara1(iefani,k),posvara2(iefani2,l))= &
                 sa(posvara1(iefani,k),posvara2(iefani2,l))+ &
                 b(pos1,k)*b(pos2,l)*dia(i) 
              enddo
             enddo
            enddo
          enddo
! elementos no diagonales
          iplace=ifirst(i)
! Salta elementos que no son A-1
1014      if((ivcol(iplace).le.ifac(nfac)).and.(iplace.gt.0))then
            iplace=inext(iplace)
            goto 1014
          endif
1015      if(iplace.gt.0)then
            do iefani=1,nefani
              pos1=i-ifac(nfac)+ifac(nfac-nefani+iefani)
              do iefani2=1,nefani
                pos2=ivcol(iplace)-ifac(nfac)+ifac(nfac-nefani+iefani2)

                do  k=1,ncar
                  do  l=1,ncar
                    sa(posvara1(iefani,k),posvara2(iefani2,l))= &
                     sa(posvara1(iefani,k),posvara2(iefani2,l))+ &
                     b(pos1,k)*b(pos2,l)*zhz(iplace)
                  enddo
                enddo
              enddo
            enddo
            iplace=inext(iplace)
            goto 1015
          endif
        enddo

        na=nanim-ngrup
        sa=sa/real(na-NCAR*nefani-1)       

!      do
        call ginv1(sa,ncar*nefani,ncar*nefani,tol,irango)
        call wish(ncar*nefani,sa,va,na-NCAR*nefani-1,x1)
        call ginv1(va,ncar*nefani,ncar*nefani,tol,irango)


!EV subrutina de denseop (Misztal)
!        call pos_def(ve(1:ncar,1:ncar),'va no positiva def, corregida',1.d-6,stat)        
!        if (stat) npdva=npdva+1

!EV lo he metido dentro de un loop. Igual para el resto de las varianzas
!		call eigen(ve(1:ncar,1:ncar),d(1:ncar),v(1:ncar,1:ncar))
!		if (minval(d(1:ncar)).gt.1.d-1) exit
!		print*, 'eigenvalueva' , minval(d(1:ncar))
!		print*, 'matriz eig', ve(1:ncar,1:ncar), d(1:ncar)
!        pause
!      enddo  
		vara=va

!	print *,'pasa7',ssumyeee(1,1,1,1)
! muestreo de los mendelian sampling
	do j=1,ncar
	   do i=ifac(nfac)+1,ifac(nfac+1)-ngrup
              do iefani=1,nefani
                pos1=i-ifac(nfac)+ifac(nfac-nefani+iefani)
                pos2=iped(i-ifac(nfac),2)+ifac(nfac-nefani+iefani)
                pos3=iped(i-ifac(nfac),3)+ifac(nfac-nefani+iefani)
!               print *,i-ifac(nfac),iped(i-ifac(nfac),2),iped(i-ifac(nfac),3),pos2,pos3
		xmend(i-ifac(nfac),j,iefani)=b(pos1,j)-0.5*b(pos2,j)-0.5*b(pos3,j)
!		print *,i-ifac(nfac),j,iefani,xmend(i-ifac(nfac),j,iefani)
              enddo
           enddo
	enddo
!	pause
!	print *,'pasa8',ssumyeee(1,1,1,1)
!	reconstrucción de breeding values por grupos
	bv=0
	sumy=0
	ssvar=0
	
	do j=1,ncar
 	   do iefani=1,nefani
 	      do i=1,nanim-ngrup
 	         do k=1,nsegm
 	         if (iped(i,4).eq.k) then
 	            bv(i,j,iefani,k)=xmend(i,j,iefani)
 	         else
 	            bv(i,j,iefani,k)=0.
 	         endif
 	         if (iped(i,2).lt.i) then
 	            bv(i,j,iefani,k)=bv(i,j,iefani,k)+0.5*bv(iped(i,2),j,iefani,k)
 	         endif
 	         if (iped(i,3).lt.i) then
 	            bv(i,j,iefani,k)=bv(i,j,iefani,k)+0.5*bv(iped(i,3),j,iefani,k)
 	         endif
 	         sumy(j,iefani,k,iped(i,5))=sumy(j,iefani,k,iped(i,5))+bv(i,j,iefani,k)
 	         enddo
 	         bv(i,j,iefani,nsegm+1)=xmend(i,j,iefani)
 	         if (iped(i,2).lt.i) then
 	            bv(i,j,iefani,nsegm+1)=bv(i,j,iefani,nsegm+1)+0.5*bv(iped(i,2),j,iefani,nsegm+1)
 	         endif
 	         if (iped(i,3).lt.i) then
 	            bv(i,j,iefani,nsegm+1)=bv(i,j,iefani,nsegm+1)+0.5*bv(iped(i,3),j,iefani,nsegm+1)
 	         endif
 	         sumy(j,iefani,nsegm+1,iped(i,5))=sumy(j,iefani,nsegm+1,iped(i,5))+bv(i,j,iefani,nsegm+1)
 	         if (ndesc(i).gt.0) then
 	           sumsel(j,iefani,iped(i,5))=sumsel(j,iefani,iped(i,5))+bv(i,j,iefani,nsegm+1)
 	           sumwei(j,iefani,iped(i,5))=sumwei(j,iefani,iped(i,5))+bv(i,j,iefani,nsegm+1)*ndesc(i)
 	         endif
 	           ssvar(j,iefani,iped(i,5))=ssvar(j,iefani,iped(i,5))+bv(i,j,iefani,nsegm+1)**2.
 	      enddo      
 	   enddo
 	enddo
	
!	print *,'pasa9',ssumyeee(1,1,1,1)
	do j=1,ncar
	do iefani=1,nefani
	do ki=1,nyear
	do k=1,nsegm
	sumy(j,iefani,k,ki)=sumy(j,iefani,k,ki)/ny(ki)
!	print *,j,iefani,k,ki,sumy(j,iefani,k,ki)
	enddo
	sumy(j,iefani,nsegm+1,ki)=sumy(j,iefani,nsegm+1,ki)/ny(ki)
!	print *,j,iefani,ki,varyear(j,iefani,ki)
!	print *,ki,sumy(j,iefani,nsegm+1,ki)
	enddo
	enddo
	enddo
	
	do j=1,ncar
	do iefani=1,nefani
	do k=1,nsegm
	write(100+j*10+iefani,'(2i10,1000f10.3)')ijk,k,(sumy(j,iefani,k,ki),ki=1,nyear)
	enddo
	k=nsegm+1
	write(100+j*10+iefani,'(i10,a10,1000f10.3)')ijk,'TOT',(sumy(j,iefani,k,ki),ki=1,nyear)
	enddo
	enddo

!	calculo de la respuesta por grupo y año
	
    endif VarianceComponentsEstimated

    if(.not.(VCE)) then
      !	put variances inverted back
      call ginv1(vare,ncar,ncar,tol,irango)
      call ginv1(vara,ncar*nefani,ncar*nefani,tol,irango)
      do i=1,nper
	se=varp(:,:,i)
	call ginv1(se,ncar,ncar,tol,irango)
	varp(:,:,i)=se
      enddo
    endif

9998  continue 

      do i=1,ncar*nefani
        print *,(varA(i,j),j=1,ncar*nefani)
      enddo

      do iper=1,nper
      do i=1,ncar
        print *,(varP(i,j,iper),j=1,ncar)
      enddo
      enddo

      do i=1,ncar
        print *,(vare(i,j),j=1,ncar)
      enddo
      print *,' imue ',ijk
      call printtime

      do i=1,ncar
        total=vara(i,i)+vare(i,i)
        do j=1,nper
          total=total+varp(i,i,j)
        enddo
        h2(i)=vara(i,i)/total
      enddo
! Var y cov. excepto cov. de permanente
!      write(15,'(20f15.8)')
!     + ((vara(i,j),j=i,ncar*nefani),i=1,ncar*nefani),
!     +((varp(i,i,j),j=1,nper),i=1,ncar),
!     +((vare(i,j),j=i,ncar),i=1,ncar),
!     +(h2(i),i=1,ncar)
! Var y cov. permanente incluido
!      write(15,'(20f15.8)')

! Write h2 if animal model only	
     if(VCE) then
       if((tipomodelo/='sire').and.(nefani==1)) then
	 write(15,'(40f15.8)')  &
	    ((vara(i,j),j=i,ncar*nefani),i=1,ncar*nefani), &
	    (((varp(i,j,k),j=i,ncar),i=1,ncar),k=1,nper),  &
	    ((vare(i,j),j=i,ncar),i=1,ncar), &
	    (h2(i),i=1,ncar)
       else	    
	 write(15,'(40f15.8)')  &
	    ((vara(i,j),j=i,ncar*nefani),i=1,ncar*nefani), &
	    (((varp(i,j,k),j=i,ncar),i=1,ncar),k=1,nper),  &
	    ((vare(i,j),j=i,ncar),i=1,ncar)
       endif
     endif
	
! ---------	
! Contrasts
! ---------
!    uncomment next line if you want contrasts
    write(20,'(20f15.8)') b(1:3,1),b(1:3,2)
! -------------
! end contrasts
! -------------

!      write(15,'(20f15.8)')((vara(i,j),j=i,ncar),i=1,ncar),
!     +(h2(i),i=1,ncar)
      write(16,'(20f15.8)')((thr(j,i),i=1,nthr(j)),j=1,ncarumb)

!EV fichero de continuaci\F3n

		       if (mod(ijk,iguard).eq.0) then
!          open(unit=14,file='GS_'//adjustl(fichero),status='replace')
             open(21, file='continuacion.txt')
				     write(21,'(a,i10)')' numero de iteracion:',ijk*icad
             write(21,'(a,i10)')' burn-in:',lap*icad
             write(21,*)
             write(21,*)' varianzas  y umbrales'
             write(21,*)  
				     write(21,'(a,i10)')'matriz (co)varianzas'
             do i=1, ncar*nefani
					     write(21,'(20f15.8)')(vara(i,j),j=1,ncar*nefani)
					   enddo
					   do k=1,nper
					     write(21, '(a,i10)') 'matriz permanentes', k
					     do i=1,ncar
						     write(21,'(20f15.8)')(varp(i,j,k),j=1,ncar)
					     enddo
					   enddo
					   write(21,'(a,i10)') 'matriz residuales'
					   do i=1,ncar
					     write(21,'(20f15.8)') (vare(i,j),j=1,ncar)
					   enddo
					   write(21,'(a,i10)') 'umbrales'
					   do j=1,ncarumb
					     write(21,'(20f15.8)') (thr(j,i),i=1,nthr(j))
					   enddo
             write(21,*) 'semilla', x1
				     write(21,*) 'b'
				     do i=1, nrow
				       write(21,*)(b(i,1:ncar))
				     enddo
             close(21)
           endif




!-- Preparacion de la salida        
      AfterBurnin: if (ijk.gt.lap) then
      
        VarianceComponentEstimation: if(VCE) then
	  do i=1,ncar
            do j=1,ncar
              sumvare(i,j)=sumvare(i,j)+vare(i,j)
              do k=1,nper
        	sumvarp(i,j,k)=sumvarp(i,j,k)+varp(i,j,k)
              enddo
              ssvare(i,j)=ssvare(i,j)+vare(i,j)**2
              do k=1,nper
        	ssvarp(i,j,k)=ssvarp(i,j,k)+varp(i,j,k)**2
              enddo
              if(i.ne.j)then
        	do k=1,nper
   !             para evitar problemas en efectos permanentes que no se usan
        	  if ((varp(i,i,k)*varp(j,j,k)).le.0.d0) then
                    corp(i,j,k)=0.d0
        	  else
                    corp(i,j,k)=varp(i,j,k)/sqrt(varp(i,i,k)*varp(j,j,k))
        	  endif
        	enddo
        	core(i,j)=vare(i,j)/sqrt(vare(i,i)*vare(j,j))
              endif
              if(i.eq.j)then
        	total=vare(i,i)
        	do k=1,nefani
                   total=total+vara(posvara1(k,i),posvara1(k,i))
        	enddo
        	do k=1,nper
        	  total=total+varp(i,i,k)
        	enddo
        	do k=1,nefani
        	  cora(posvara1(k,i),posvara1(k,i))= &
        	  vara(posvara1(k,i),posvara1(k,i))/total
        	enddo
        	do k=1,nper
        	  corp(i,j,k)=varp(i,i,k)/total
        	enddo
        	core(i,j)=vare(i,i)/total
              endif
              do k=1,nper
        	sumcorp(i,j,k)=sumcorp(i,j,k)+corp(i,j,k)
              enddo
              sumcore(i,j)=sumcore(i,j)+core(i,j)
              do k=1,nper
        	sscorp(i,j,k)=sscorp(i,j,k)+corp(i,j,k)**2
              enddo
              sscore(i,j)=sscore(i,j)+core(i,j)**2
            enddo
	  enddo
	  do i=1,ncar*nefani
            do j=1,ncar*nefani
              sumvara(i,j)=sumvara(i,j)+vara(i,j)
              ssvara(i,j)=ssvara(i,j)+vara(i,j)**2
              if(i.ne.j)then
        	if ( (vara(i,i)*vara(j,j)).le.0.d0 )then
        	  cora=0.d0
        	else  
        	  cora(i,j)=vara(i,j)/sqrt(vara(i,i)*vara(j,j))
        	endif
              endif
              sumcora(i,j)=sumcora(i,j)+cora(i,j)
              sscora(i,j)=sscora(i,j)+cora(i,j)**2
            enddo
	  enddo
        endif VarianceComponentEstimation
	

	
	!mean and sd of solutions
	do i=1,nrow
	  do j=1,ncar
	    sumb(i,j)=sumb(i,j)+b(i,j)
	    ssb(i,j)=ssb(i,j)+b(i,j)**2
	  enddo
	enddo
	
	!mean and sd of genetic trends
	
!	ssumyeee=0
	
	do j=1,ncar
	do iefani=1,nefani
	do k=1,nsegm+1
	do ki=1,nyear
!	print *,j,iefani,k,ki,sumy(j,iefani,k,ki),ssumyeee(j,iefani,k,ki)
	ssumyeee(j,iefani,k,ki)=ssumyeee(j,iefani,k,ki)+sumy(j,iefani,k,ki)
	sssumyeee(j,iefani,k,ki)=sssumyeee(j,iefani,k,ki)+(sumy(j,iefani,k,ki))**2.
	print *,ssumyeee(j,iefani,k,ki),sssumyeee(j,iefani,k,ki),sumy(j,iefani,k,ki),ki,k,iefani,j
	enddo
	enddo
	enddo
	enddo
	

        SavingResults:  if (mod(ijk-lap,100).eq.0) then




  !          open(unit=14,file='GS_'//adjustl(fichero),status='replace')
          if(VCE) then

!	  escribe las tendencias genéticas


	   
	    open(173,file='trends.txt',status='replace')
	    do j=1,ncar
   	       do iefani=1,nefani
	         do k=1,nsegm+1
	           do ki=1,nyear
                     xmean=ssumyeee(j,iefani,k,ki)/(ijk-lap)
                     sd=sssumyeee(j,iefani,k,ki)-ssumyeee(j,iefani,k,ki)*ssumyeee(j,iefani,k,ki)/(ijk-lap)
	             sd=sqrt(sd/(ijk-lap-1))
	             write(173,*)j,iefani,k,ki,xmean,sd,ny(ki)
		   enddo
	          enddo
	        enddo
	     enddo
	     close(173)
	     
!	escribe la evalución de las varianzas genéticas

    	  
            open(unit=14,file='results.txt',status='replace')
            write(14,'(a,a)')'Parameter file: ',fichero
            write(14,'(a,i10)')' Iteration number:',ijk*icad
            write(14,'(a,i10)')' Burn-in:',lap*icad
            write(14,*)
		    write(14,*)'ve stat=true', npdve
		    write(14,*)'vp stat=true', npdvp
		    write(14,*)'va stat=true', npdva
            write(14,*)'Average additive variance'
            do i=1,ncar*nefani
              write(14,'(20f15.8)')((sumvara(i,j)/(ijk-lap)), &
                     j=1,ncar*nefani)
            enddo
            write(14,*)' Sd Additive variance'
            do i=1,ncar*nefani           
              write(14,'(20f15.8)') (sqrt1((ssvara(i,j)-(sumvara(i,j)**2)/ &
               (ijk-lap)) / (ijk-lap-1)),j=1,ncar*nefani)
            enddo
            do k=1,nper
              write(14,*)
              write(14,'(a,i3,a)')' Average environmental variance  ',k,'-th'
              do i=1,ncar
        	write(14,'(20f15.8)')((sumvarp(i,j,k)/(ijk-lap)),j=1,ncar)
              enddo
              write(14,*)' Sd environmental variance'
              do i=1,ncar
        	write(14,'(20f15.8)') (sqrt1((ssvarp(i,j,k)- &
        	 (sumvarp(i,j,k)**2)/ &
               (ijk-lap)) / (ijk-lap-1)),j=1,ncar)
              enddo
            enddo
            write(14,*)
            write(14,*)' Average residual variance'
            do i=1,ncar
              write(14,'(20f15.8)')((sumvare(i,j)/(ijk-lap)),j=1,ncar)
            enddo
            write(14,*)' Sd residual variance '
            do i=1,ncar
              write(14,'(20f15.8)') (sqrt1((ssvare(i,j)-(sumvare(i,j)**2)/ &
             (ijk-lap)) / (ijk-lap-1)),j=1,ncar)
            enddo
	    
	    NOTWEIRDMODEL: if((tipomodelo/='sire').and.(nefani==1)) then
	    
            write(14,*)
            write(14,*)' Average h2 and additive correlation '
            do i=1,ncar*nefani
              write(14,'(20f15.8)')((sumcora(i,j)/(ijk-lap)), &
                  j=1,ncar*nefani)
            enddo
            write(14,*)' Sd h2 and additive correlation'
            do i=1,ncar*nefani
              write(14,'(20f15.8)') (sqrt1((sscora(i,j)-(sumcora(i,j)**2)/ &
              (ijk-lap)) / (ijk-lap-1)),j=1,ncar*nefani)
            enddo
            do k=1,nper
              write(14,*)
              write(14,'(a,i3,a)')' Average c2 and environmental cor ',k,'-th'
              do i=1,ncar
        	write(14,'(20f15.8)')((sumcorp(i,j,k)/(ijk-lap)),j=1,ncar)
              enddo
              write(14,*)' Sd c2 and environmental cor'
              do i=1,ncar
        	write(14,'(20f15.8)') (sqrt1((sscorp(i,j,k)- &
        	 (sumcorp(i,j,k)**2)/ &
        	 (ijk-lap)) / (ijk-lap-1)),j=1,ncar)
              enddo
            enddo
            write(14,*)
            write(14,*)' Average he2 and residual cor'
            do i=1,ncar
              write(14,'(20f15.8)')((sumcore(i,j)/(ijk-lap)),j=1,ncar)
            enddo
            write(14,*)' Sd he2 and residual cor'
            do i=1,ncar
              write(14,'(20f15.8)') (sqrt1((sscore(i,j)-(sumcore(i,j)**2)/ &
               (ijk-lap)) / (ijk-lap-1)),j=1,ncar)
            enddo
	    ELSE
	    
	      WRITE(14,*) '--------------------- NOTE -----------------------------------'
	      WRITE(14,*) 'Heritabilities and r_g have to be calculated from ''samples.txt'' '
	      WRITE(14,*) '--------------------------------------------------------------'
	      
	    ENDIF NOTWEIRDMODEL

            close(14)
	  endif


          open(17,file='solutions.txt',status='replace')
          do i=1,nrow
              write(17,'(20f15.8)')(sumb(i,j)/(ijk-lap), &
     sqrt1( (ssb(i,j)-(sumb(i,j)**2)/ &
            (ijk-lap)) / (ijk-lap-1) ),j=1,ncar)
          enddo
          close(17)


        endif SavingResults


      endif Afterburnin

9999  continue 
!-- Fin del gibbs sampler

call print_version(start=.false.)

      stop 
      end 

