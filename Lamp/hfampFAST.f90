!****************************************************
!
! Angular momentum projection of slater determinants
!
!    Dr. Calvin Johnson
!    Joshua Staker
!
!****************************************************

! Main integration routine
!
!    alp, bet, gam - indicies for quadrature integration
!    js - User defined angular momentum (coupled - integer (or half-integer) value)
!    np - dimension of generated H and N matrixes (np = 2*js + 1)
!    mm, mp - indices for generated matrix
!    hmm, hmp - actual values of m', m for wignerD matrix generation
!    uw - Generated WignerD matrix from user js
!
!    jerkap, jerkbp, jerkgp - values to store integration
!        producing H matrix (Dimension np X np)
!
!    jerkapn, jerkbpn, jerkgpn - values to store integration
!        producing N (normalization) matrix (Dimension np X np)
!
!    psd, nsd - proton/neutron slater determinants
!    psdr, nsdr - Rotated slater determinants via RotMat matrix
!    RotMat - Block diagonal rotation matrix generated in
!        subroutine Psi_New
!


SUBROUTINE JackedTest(js, np, psd, nsd, nf, pevals,posdef,ntrace,htrace)

USE phf_vals
USE system_parameters
USE spstate
USE psis

IMPLICIT NONE

INTEGER (KIND = 8), INTENT(IN) :: np
REAL (KIND = 8), INTENT(IN) :: js
COMPLEX (KIND = 8), INTENT(IN) :: psd(nsps,numprot), nsd(nsps,numneut)
LOGICAL, INTENT(INOUT) :: posdef
INTEGER (KIND = 8), INTENT(OUT) :: nf
REAL (KIND = 8), INTENT(OUT) :: ntrace,htrace
REAL (KIND = 8), INTENT(OUT) :: pevals(np)

INTEGER :: it, ii,jj,kk,ctest
INTEGER (KIND = 8) :: mp,mm,npts,nptsl,alp,bet,gam,errorflag,loop1,loop2,tt,hcheck
INTEGER, ALLOCATABLE, DIMENSION(:) :: skip
REAL :: tol
REAL (KIND = 8) :: i,a,b,hmp,hm,k,lam,pi,ehartree,ooeps,tol_test = 1.0d-5
REAL (KIND = 8) :: hbool
REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: xleg,wleg,xlegb,wlegb
REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:,:) :: xjac,wjac
COMPLEX (KIND = 8) :: uw!,jerkap,jerkapn,jerkbp,jerkbpn,jerkgp,jerkgpn
COMPLEX (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: phfA,phfB,phfG,NphfA,NphfB,NphfG,IntWig
COMPLEX (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: RotMat, L, Linv,RotTesta,RotTestb,RotTestg
LOGICAL :: invFlag = .true.
PARAMETER(pi = 2.0d0*DASIN(1.0d0))
PARAMETER(ooeps = 1.0d0/(8.0d0*pi*pi))

COMPLEX (KIND = 8),ALLOCATABLE :: psdr(:,:),nsdr(:,:)
COMPLEX (KIND = 8),ALLOCATABLE :: rhopij(:,:),rhonij(:,:),Hmat(:,:),Nmat(:,:),Hprime(:,:)
COMPLEX (KIND = 8),ALLOCATABLE :: Hmatp(:,:), Hmatn(:,:), Nmatp(:,:), Nmatn(:,:)
COMPLEX (KIND = 8) :: ovlpp,ovlpn
COMPLEX (KIND = 8) :: vme,ic

integer :: info,lwork
real (kind = 8), ALLOCATABLE :: e(:)
real (kind = 8),ALLOCATABLE :: rwork(:)
COMPLEX (kind = 8), ALLOCATABLE :: work(:)
REAL (KIND = 8), ALLOCATABLE :: pairOP(:)

1001 format(A,1x,f7.3,1x,A,1x,g15.8)

ALLOCATE (psdr(nsps,numprot),nsdr(nsps,numneut))
ALLOCATE (pairOP(nsps))

!CALL GetSherpaSD(psd,nsd)

hcheck = 0
i = 2
DO
    hcheck = hcheck + 1
    i = i + 2*hcheck
    IF (i >= nsps) THEN
        EXIT
    END IF
END DO
    
k = 1
DO ii = 1, hcheck
    tt = spsqn(kk)%j
    kk = kk+tt+1
END DO

!write(*,*) "Enter grid dim:"
!read(*,*) ctest
ctest = 20
npts = ctest
nptsl = ctest
hbool = (4.0d0*DASIN(1.0d0))/DBLE(nptsl)
tol = 0.1e0
ic = CMPLX(0.0d0,1.0d0)

ALLOCATE(xleg(nptsl),wleg(nptsl))

ALLOCATE(xlegb(npts),wlegb(npts))
ALLOCATE(xjac(np,np,npts),wjac(np,np,npts))

ALLOCATE(RotTesta(nsps,nsps),RotTestb(nsps,nsps),RotTestg(nsps,nsps))

ALLOCATE(Hmatp(np,np),Hmatn(np,np),Nmatp(np,np),Nmatn(np,np))
ALLOCATE(Hmat(np,np),Nmat(np,np),Hprime(np,np))
ALLOCATE(phfA(np,np),phfB(np,np),phfG(np,np),NphfA(np,np),NphfB(np,np),NphfG(np,np))
ALLOCATE(IntWig(np,np))

ALLOCATE(RotMat(nsps,nsps))
ALLOCATE (rhopij(nsps,nsps),rhonij(nsps,nsps))
!Compute weights and abscissas for the alpha and gamma
!integrals using Gauss-Legendre quadrature

CALL gauleg(0.0d0,2.0d0*DASIN(1.0D0),xlegb,wlegb,nptsl)
CALL gauleg(0.0d0,4.0d0*DASIN(1.0d0),xleg,wleg,nptsl)
!CALL Trap_Points(0.0d0,4.0d0*DASIN(1.0d0),nptsl,xleg,wleg)

!NOTE: If using boole's rule, be sure to change alp and gam starting indicies to 0
!CALL bool_points(nptsl,hbool,0.0d0,4.0d0,wleg,xleg)


Hmat = CMPLX(0.0d0,0.0d0)
Nmat = CMPLX(0.0d0,0.0d0)
Hprime = CMPLX(0.0d0,0.0d0)


!!$OMP PARALLEL DO
        !Main integration routine, integrates first by alpha, then beta, then gamma

        phfA = CMPLX(0.0d0,0.0d0)
        NphfA = CMPLX(0.0d0,0.0d0)
        !RotTesta = CMPLX(0.0d0,0.0d0)
        DO alp = 1, nptsl
            phfB = CMPLX(0.0d0,0.0d0)
            NphfB = CMPLX(0.0d0,0.0d0)
            !RotTestb = CMPLX(0.0d0,0.0d0)
            DO bet = 1, nptsl
                phfG = CMPLX(0.0d0,0.0d0)
                NphfG = CMPLX(0.0d0,0.0d0)
                !RotTestg = CMPLX(0.0d0,0.0d0)
                DO gam = 1, nptsl

                    !  NOTE:  CURRENTLY SET FOR GAUSS-LEGENDRE QUADRATURE (SQUARE GRID) 

		          CALL Wigner_d2(js,int(np,kind=4),xleg(alp),xlegb(bet),xleg(gam),IntWig)
                    ! Generate block diagonal matrix, use it to rotate slater dets.
                    ! and generate rho function
                    CALL Psi_New(xleg(alp),xlegb(bet),xleg(gam),RotMat)
		    
                
                    !Using generated R matrix, Rotate the neutron and proton slater determinants
                    
				psdr = MATMUL(RotMat,psd)
                    nsdr = MATMUL(RotMat,nsd)

                    ! Generate neutron and proton overlaps (ovlp*) as well as the density matrices (roh*ij)

                    CALL makerhoij(1,numprot,psdr,psd,ovlpp,rhopij)
                    CALL makerhoij(2,numneut,nsdr,nsd,ovlpn,rhonij)

                    !Generate <H> = vme for each step in the integration

                    CALL TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
				    !write(*,'(1x,100G14.5)') xleg(alp),hm,hmp,ovlpp,ovlpn

		    		phfG = phfG + wleg(gam)*IntWig*vme*DSIN(xlegb(bet))
				    NphfG  = NphfG + wleg(gam)*IntWig*ovlpp*ovlpn*DSIN(xlegb(bet))
                    !jerkgp = jerkgp + wleg(gam)*EXP(ic*hm*xleg(alp))*(uw)*EXP(ic*hmp*xleg(gam))*vme
                    !jerkgpn = jerkgpn + wleg(gam)*EXP(ic*hm*xleg(alp))*(uw)*EXP(ic*hmp*xleg(gam))*ovlpp*ovlpn                    
                    !Test the rotation Matrix R via integration of |R|^2 using Gauss-Legendre
                    !CALL Psi_New(xleg(alp),xlegb(bet),xleg(gam),RotMat)
                    !RotTestg = RotTestg + wleg(gam)*MATMUL(DCONJG(RotMat),RotMat)*DSIN(xlegb(bet))
                END DO
                phfB = phfB + wlegb(bet)*phfG
                NphfB = NphfB + wlegb(bet)*NphfG
                !jerkbp = jerkbp + wjac(mp,mm,bet)*jerkgp
                !jerkbpn = jerkbpn + wjac(mp,mm,bet)*jerkgpn
                !RotTestb = RotTestb + wlegb(bet)*RotTestg
            END DO
			phfA = phfA + wleg(alp)*phfB
			NphfA = NphfA + wleg(alp)*NphfB
            !jerkap = jerkap + wleg(alp)*jerkbp
            !jerkapn = jerkapn + wleg(alp)*jerkbpn
            !RotTesta = RotTesta + wleg(alp)*RotTestb
        END DO
        Hmatp = phfA
		!Hmatp(mm,mp) = CONJG(Hmatp(mp,mm))
        Nmatp = NphfA
		!Nmatp(mm,mp) = CONJG(Nmatp(mp,mm))
!!$OMP END DO

Hmat = Hmatp
phf(jloc)%h = Hmat

Nmat = Nmatp
phf(jloc)%n = Nmat
!
!	Hermitian test
!
!do ii = 1, np
	!if(dble(Nmat(ii,jj))-dble(Nmat(jj,ii)) < tol_test) write(*,*) "You win"
!	write(*,'(1x,100G14.5)') (dble(Nmat(ii,jj)),jj=1,np)
!end do
!
!	J contribution Test
!
ntrace = 0.0d0
htrace = 0.0d0
do it = 1, np
    ntrace = ntrace + dble(Nmat(it,it))
    htrace = htrace + dble(Hmat(it,it))
end do
ntrace = ntrace*dble(np)*ooeps
write(*,1001) "J =",js,"contribution:",ntrace

nf = np
ALLOCATE(L(np,np),Linv(np,np),skip(np))
CALL choleskymaster(np,np,Nmat,tol,invFlag,L,Linv,skip,posdef)

IF (posdef .EQV. .true.) THEN
    CALL choleskyinvert(np,np,nf,Hmat,Hprime,Linv,skip) ! Returns H N^-1
    !WRITE(*,*) nf
ELSE
    WRITE(*,*) 'Matrix not positive-definite, skipping'
    RETURN
ENDIF
! If no states are kept by Cholesky, then move on to next J state
IF (nf == 0) THEN
    RETURN
ENDIF



ALLOCATE(e(nf))
lwork = 2*nf - 1
ALLOCATE(work(lwork))
ALLOCATE(rwork(3*nf -2  ))
!WRITE(*,*) nf
CALL zheev('n','u',nf,Hprime,np,e,work,lwork,rwork,info)

DO loop1 = 1, nf
    pevals(loop1) = e(loop1)
END DO

END SUBROUTINE JackedTest

SUBROUTINE UserWigner(a,b,k,lam,js,bet,uw)

IMPLICIT NONE

REAL (KIND = 8), INTENT(IN) :: a,b,k,lam,js,bet
COMPLEX (KIND = 8), INTENT(OUT) :: uw

REAL (KIND = 8) :: kpbb,jmkkpa,jacobi,sctol,sc,st,ct
PARAMETER(sctol = 1.0d-7)

CALL ChooseFac(k+b,b,kpbb)
CALL ChooseFac(2.0d0*js-k,k+a,jmkkpa)
CALL JacobiP(INT(k),INT(a),INT(b),bet,jacobi)

uw = ((-1.0d0)**lam)*(0.5d0**((INT(a+b))/2.0d0))*DSQRT(jmkkpa/kpbb)*jacobi

END SUBROUTINE UserWigner

SUBROUTINE Factorial(n,fact)

IMPLICIT NONE

REAL (KIND = 8), INTENT(IN) :: n        !Integer to be factorialized
REAL (KIND = 8), INTENT(OUT) :: fact        !Variable to be exported -- factorial = n!
REAL (KIND = 8) :: temp, d_i            ! i -- index for DO loop, temp -- dummy variable for accumulation of n!
INTEGER :: i,i_n

temp = DLOG(DBLE(1.0))

i_n = INT(n)
DO i = 1,i_n
    d_i = DBLE(i)
    temp = temp+DLOG(d_i)
END DO

fact = temp

END SUBROUTINE Factorial

!Choose Fac generates binomial coefficients
!    i.e. (a choose b) = a!/(b!(a - b)!)

SUBROUTINE ChooseFac(a,b,fact)

IMPLICIT NONE

REAL (KIND = 8), INTENT(IN) :: a,b
REAL (KIND = 8), INTENT(OUT) :: fact

REAL (KIND = 8) :: afact,bfact,ambfact

CALL Factorial(a,afact)
CALL Factorial(b,bfact)
CALL Factorial(a-b,ambfact)

fact = DEXP(afact - bfact - ambfact)

END SUBROUTINE ChooseFac

!Jacobi Polynomial subroutine

SUBROUTINE JacobiP(n,alp,bet,x,Jac)

IMPLICIT NONE

INTEGER, INTENT(IN) :: n, alp, bet
REAL (KIND = 8), INTENT(IN) :: x
REAL (KIND = 8), INTENT(OUT) :: Jac

INTEGER :: i,j
REAL (KIND = 8) :: npalp,npbet,nu,npalpmnu,nmnu,betpnu
REAL (KIND = 8) :: temp,di

CALL Factorial(DBLE(n+alp),npalp)
CALL Factorial(DBLE(n+bet),npbet)

temp = 0
DO i = 0,n
    di = DBLE(i)
    CALL Factorial(di,nu)
    CALL Factorial(DBLE(n+alp-i),npalpmnu)
    CALL Factorial(DBLE(n-i),nmnu)
    CALL Factorial(DBLE(bet+i),betpnu)
    temp = temp + (1.0d0/DEXP(nu + npalpmnu + nmnu + betpnu))*((x-1.0d0)**(n-i))*(x+1.0d0)**i
END DO

Jac = (2.0d0**(-n))*DEXP(npalp+npbet)*temp

END SUBROUTINE JacobiP

!Subroutine for generating the block diagonal matrix R (called RotMat)
!
!    INPUT:
!    alpha, beta, gamma: euler angle values from integration routine
!
!    OUTPUT:
!    RotMat: Double Complex rotation matrix of dimension nsps X nsps
!
!    NOTE: as in the main integration subroutine, js corresponds
!        to the angular momentum number, however in this routine
!        js comes from the slater determinant for the fermions and
!        is therefore a half integer.

SUBROUTINE Psi_New(alpha,beta,gamma,RotMat)

USE system_parameters
USE spstate

IMPLICIT NONE

REAL (KIND = 8), INTENT(IN) :: alpha,beta,gamma
COMPLEX (KIND = 8), DIMENSION(nsps,nsps), INTENT(OUT) :: RotMat

INTEGER :: i,k,mm,mp,temp, temp1, count, g, h, ndim,tt
REAL (KIND = 8) :: js
COMPLEX (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: wignerD
 
temp = 0

k = 1
i = spsqn(k)%j+1
DO
    tt = spsqn(k)%j
    !write(*,*) k,spsqn(k)%j
    ndim = tt + 1
    js = DBLE(tt)/2.0d0
    ALLOCATE(wignerD(1:ndim,1:ndim))
    CALL Wigner_d2(js,ndim,alpha,beta,gamma,wignerD)
    DO g = 1, ndim
        temp = temp + 1
        temp1 = 0
        DO h = 1,nsps
            IF ((h > i-ndim).AND.(h <= i)) THEN
                temp1 = temp1 + 1
                RotMat(temp,h) = wignerD(g,temp1)
            ELSE
                RotMat(temp,h) = CMPLX(0.0d0,0.0d0)
            END IF
        END DO
    END DO
    DEALLOCATE(wignerD)
    k = k + tt + 1
    IF (k > nsps) THEN
        EXIT
    END IF
    i = i + spsqn(k)%j + 1
    IF (i > nsps) THEN
        EXIT
    END IF
END DO

END SUBROUTINE Psi_New

!Subroutine that uses the wigner_d function for little d
!    and adds in the exponential terms for Jz to create
!    the big D matrix

SUBROUTINE Wigner_d2(js,np,alpha,beta,gamma,wignerD)

IMPLICIT NONE

INTEGER, INTENT(IN) :: np
REAL (KIND = 8), INTENT(IN) :: js, alpha,beta,gamma
COMPLEX (KIND = 8), DIMENSION(1:np,1:np), INTENT(OUT) :: wignerD

INTEGER :: mp,mm,lam
REAL (KIND = 8) :: k,hm,hmp,jmkkpa,kpbb,a,b
REAL (KIND = 8) :: jacobi
COMPLEX (KIND = 8) :: alphac,gammac,ic, hmpc, hmc
COMPLEX (KIND = 8), DIMENSION(1:np, 1:np) :: alpexp,gamexp,smalld

ic = CMPLX(0.0d0,1.0d0)
alphac = CMPLX(alpha,0.0d0)
gammac = CMPLX(gamma,0.0d0)

!hmp, hm correspond to the m', m values of the wigner D matrix
!and range from -js, js
!
!mp, m and the integer values corresponding to the array
!indicies when the matrix is constructed where
!np is the dimension of the wigner matrix
!with a value of np = 2*js + 1

hmp = js
DO mp = 1, np
    hmpc = CMPLX(hmp,0.0d0)
    hm = js
    DO mm = 1, np
        hmc = CMPLX(hm,0.0d0)
        k = min(js - hm, js+hm, js-hmp, js+hmp)
        IF (k == js+hm) THEN
            a = hmp - hm
            lam = hmp - hm
        ELSE IF (k == js-hm) THEN
            a = hm - hmp
            lam = 0.0d0
        ELSE IF (k == js + hmp) THEN
            a = hm - hmp
            lam = 0.0d0
        ELSE
            a = hmp - hm
            lam = hmp - hm
        END IF
        b = 2.0d0*js - 2.0d0*k - a
        CALL ChooseFac(2.0d0*js - k, k + a,jmkkpa)
        CALL ChooseFac(k+b,b,kpbb)
        CALL JacobiP(INT(k),INT(a),INT(b),DCOS(beta),jacobi)
        smalld(mp,mm) = ((-1)**lam)*DSQRT(jmkkpa/kpbb)*(DSIN(beta/2.0d0)**a)*(DCOS(beta/2.0d0)**b)*jacobi
        IF (hmp == hm) THEN
            alpexp(mp,mm) = EXP(-ic*alpha*hm)
            gamexp(mp,mm) = EXP(-ic*gamma*hmp)
        ELSE
            alpexp(mp,mm) = CMPLX(0.0d0,0.0d0)
            gamexp(mp,mm) = CMPLX(0.0d0,0.0d0)
        END IF
        !Reduce m' value by 1
        hm = hm - 1.0d0
    END DO
    !Reduce m value by 1
    hmp = hmp - 1.0d0
END DO

wignerD = MATMUL(alpexp,MATMUL(smalld,gamexp))

END SUBROUTINE Wigner_d2 

SUBROUTINE Trap_Points(xstart, xend, nptstrap, xtrap, wtrap)

IMPLICIT NONE

INTEGER (KIND = 8), INTENT(IN) :: nptstrap
REAL (KIND = 8), INTENT(IN) :: xstart, xend
REAL (KIND = 8) :: dx
REAL (KIND = 8), DIMENSION(nptstrap), INTENT(OUT) :: xtrap,wtrap

INTEGER (KIND = 8) :: i

dx = (xend - xstart) / DBLE(nptstrap - 1)

wtrap = 1.0d0
wtrap(1) = 0.5d0
wtrap(nptstrap) = 0.5d0
wtrap = dx*wtrap
xtrap(1) = xstart
xtrap(nptstrap) = xend
!write(*,*) wtrap(1),xtrap(1)
DO i = 2, nptstrap-1
    xtrap(i) = xstart+DBLE(i-1)*dx
    !write(*,*) wtrap(i),xtrap(i)
END DO
!write(*,*) wtrap(nptstrap),xtrap(nptstrap)
!STOP

END SUBROUTINE Trap_Points

subroutine bool_points(N,h,a,b,w,x)

implicit none

integer, intent(in) :: N
real (kind = 8), intent(in) :: a,b,h
real (kind = 8), dimension(0:N), intent(out) :: w,x

integer :: i

x(0) = a
x(N) = b
w(0) = 7.0d0
w(N) = 7.0d0

do i = 1,N-1
	x(i) = a + i*h
	if ((mod(i,4) == 1).OR.(mod(i,4) == 3)) then
		w(i) = 32.0d0
	else if (mod(i,4) == 2) then
		w(i) = 12.0d0
	else
		w(i) = 14.0d0
	end if
end do

w = (2.0d0/45.0d0)*h*w

end subroutine bool_points
