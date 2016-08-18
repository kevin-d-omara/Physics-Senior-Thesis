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
!  subroutines called:
!    MakeParityOp
!    gauleg
!    Wigner_d2
!    Psi_new
!    makerhoij
!    TBMEmaster
!    choleskymaster
!    choleskyinvert
!    zheev
!


SUBROUTINE Projection_with_Parity(js,np,psd,nsd,tol,ParityTest,alpha_i,beta_j, &
        gamma_k,N_ijk,H_ijk,PN_ijk,PH_ijk,allOvlpp,allOvlpn,isNorm)

USE phf_vals
USE system_parameters
USE spstate
USE psis
USE errortests

IMPLICIT NONE

interface
    SUBROUTINE Psi_New(alpha,beta,gamma,RotMat)  ! INTERFACE
        USE system_parameters
        USE spstate
        IMPLICIT NONE
        REAL (KIND = 8), INTENT(IN) :: alpha,beta,gamma
        COMPLEX (KIND = 8), DIMENSION(nsps,nsps), INTENT(OUT) :: RotMat
    end subroutine Psi_New  ! INTERFACE


    SUBROUTINE Wigner_d2(js,np,alpha,beta,gamma,wignerD)  ! INTERFACE
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: np
        REAL (KIND = 8), INTENT(IN) :: js, alpha,beta,gamma
        COMPLEX (KIND = 8), DIMENSION(1:np,1:np), INTENT(OUT) :: wignerD
    end subroutine Wigner_d2    ! INTERFACE

    subroutine makerhoij(it,np,sdf,sdi,ovlp,rhoij)  ! INTERFACE
        use spstate
        use system_parameters
        implicit none
        integer it
        integer np
        complex(kind = 8) :: sdf(nsps,np),sdi(nsps,np)
        complex(kind = 8) :: rhoij(nsps,nsps)
        complex(kind = 8) :: ovlp
    end subroutine makerhoij  ! INTERFACE

    subroutine geneigsolverdiag(ndim,nparity,ipar,NMat,HMat,tol,posdef,nf,pevals,normvals)  ! INTERFACE
        use errortests
        implicit none
        integer (kind=4) :: ndim, nparity   ! dimension of matrices
        integer ipar
        COMPLEX (KIND = 8):: Hmat(ndim,ndim),NMat(ndim,ndim)   !Input hamiltonian and norm matrices
        LOGICAL, INTENT(INOUT) :: posdef(2)
        REAL :: tol
        INTEGER (KIND = 8), INTENT(OUT) :: nf
        REAL (KIND = 8), INTENT(OUT) :: pevals(nparity,ndim)
        real (kind=8) :: normvals(nparity,ndim)
    end subroutine geneigsolverdiag    ! INTERFACE
end interface

INTEGER (KIND = 8), INTENT(IN) :: np
REAL (KIND = 8), INTENT(IN) :: js
COMPLEX (KIND = 8), INTENT(IN) :: psd(numsd,nsps,numprot), nsd(numsd,nsps,numneut)
LOGICAL, INTENT(IN) :: ParityTest
LOGICAL :: posdef(2)
INTEGER (KIND = 8) :: nf
REAL (KIND = 8) :: pevals(npair,numsd*np)
real (kind=8) :: normvals(npair,numsd*np)
real (kind=8) :: ntrace(2)

integer (kind=4) :: nftmp
INTEGER (KIND = 8) :: sditi,sditj
INTEGER (KIND = 8) :: mp,mm,npts,nptsl,alp,bet,gam,errorflag,loop1,loop2,tt,hcheck,nkeep,itt
INTEGER, ALLOCATABLE, DIMENSION(:) :: skip
REAL :: tol
REAL (KIND = 8) :: a,b,hmp,hm,lam,pi,dx
REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: xleg,wleg,xlegb,wlegb
REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:,:) :: xjac,wjac
COMPLEX (KIND = 8) :: uw,jerkap,jerkapn,jerkbp,jerkbpn,jerkgp,jerkgpn,ic
COMPLEX (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: phfA, phfB, phfG, NphfA, NphfB, NphfG
COMPLEX (KIND = 8) :: Pjerkap,Pjerkapn,Pjerkbp,Pjerkbpn,Pjerkgp,Pjerkgpn
COMPLEX (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: PphfA, PphfB, PphfG, PNphfA, PNphfB, PNphfG
COMPLEX (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: IntWig
COMPLEX (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: RotMat, L, Linv,RotTesta,RotTestb,RotTestg
LOGICAL :: invFlag = .true.
PARAMETER(pi = 2.0d0*DASIN(1.0d0))
integer :: i,it,k,jtt

COMPLEX (KIND = 8),ALLOCATABLE :: psdr(:,:),nsdr(:,:),psdpr(:,:),nsdpr(:,:),ExpMat(:,:,:)
COMPLEX (KIND = 8),ALLOCATABLE :: rhopij(:,:),rhonij(:,:),Hmat(:,:),Nmat(:,:),Hprime(:,:)
COMPLEX (KIND = 8),ALLOCATABLE :: Hmatp(:,:), Hmatn(:,:), Nmatp(:,:), Nmatn(:,:)
COMPLEX (KIND = 8),ALLOCATABLE :: PHmatp(:,:), PHmatn(:,:), PNmatp(:,:), PNmatn(:,:)
!Multiple Slater Determinant matrices
COMPLEX (KIND = 8),ALLOCATABLE :: SDHmat(:,:),SDNmat(:,:),SDPHmat(:,:),SDPNmat(:,:)
COMPLEX (KIND = 8),ALLOCATABLE :: TSDHmat(:,:),TSDNmat(:,:)
COMPLEX (KIND = 8) :: ovlpp,ovlpn
COMPLEX (KIND = 8) :: vme

integer :: info,lwork,pairP,pairOUT,ii
logical :: testing
real (kind = 8) :: trnorm
real (kind = 8), ALLOCATABLE :: e(:)
real (kind = 8),ALLOCATABLE :: rwork(:)
COMPLEX (kind = 8), ALLOCATABLE :: work(:)
REAL (KIND = 8), ALLOCATABLE :: parOP(:), pairMat(:,:)
real (kind=8) :: projfactor   ! (2j+1)/8pi^2, in order to correctly project out

real (kind=8), allocatable, intent(in) :: alpha_i(:), gamma_k(:), beta_j(:)
complex (kind = 8), allocatable, intent(inout) :: N_ijk(:,:,:), H_ijk(:,:,:), PN_ijk(:,:,:), PH_ijk(:,:,:)
complex (kind = 8), allocatable, intent(inout) :: allOvlpp(:,:,:,:), allOvlpn(:,:,:,:)
logical, intent(in) :: isNorm

npts = int(2.d0*js)+1 !2*Jmax+1

ALLOCATE (psdr(nsps,numprot),nsdr(nsps,numneut))
ALLOCATE (parOP(nsps),pairMat(nsps,nsps))

ic = CMPLX(0.0d0,1.0d0)
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
DO i = 1, hcheck
    tt = spsqn(k)%j
    !WRITE(*,*) tt/2.0d0
    k = k+tt+1
END DO


IF (ParityTest) CALL MakeParityOp(parOP)

!npts = 35  ! set from outside
!tol = 0.1e0  ! set from outside
ALLOCATE(xleg(npts),wleg(npts))
ALLOCATE(xlegb(npts),wlegb(npts))
ALLOCATE(xjac(np,np,npts),wjac(np,np,npts))

ALLOCATE(RotTesta(nsps,nsps),RotTestb(nsps,nsps),RotTestg(nsps,nsps))

ALLOCATE(Hmatp(np,np),Hmatn(np,np),Nmatp(np,np),Nmatn(np,np))
ALLOCATE(PHmatp(np,np),PHmatn(np,np),PNmatp(np,np),PNmatn(np,np))
ALLOCATE(Hmat(np,np),Nmat(np,np),Hprime(numsd*np,numsd*np))

ALLOCATE(RotMat(nsps,nsps),ExpMat(nsps,nsps,npts ))
ALLOCATE (rhopij(nsps,nsps),rhonij(nsps,nsps))

ALLOCATE(phfA(np,np),phfB(np,np),phfG(np,np),NphfA(np,np),NphfB(np,np),NphfG(np,np))
ALLOCATE(PphfA(np,np),PphfB(np,np),PphfG(np,np),PNphfA(np,np),PNphfB(np,np),PNphfG(np,np))
ALLOCATE(IntWig(np,np))

ALLOCATE(SDHmat(numsd*np,numsd*np),SDNmat(numsd*np,numsd*np),SDPHmat(numsd*np,numsd*np),SDPNmat(numsd*np,numsd*np))
ALLOCATE(TSDHmat(numsd*np,numsd*np),TSDNmat(numsd*np,numsd*np))

!Compute weights and abscissas for the alpha and gamma
!integrals using Gauss-Legendre quadrature

CALL gauleg(0.0d0,4.0d0*DASIN(1.0d0),xleg,wleg,npts)
CALL gauleg(0.0d0,2.0d0*DASIN(1.0d0),xlegb,wlegb,npts)
!CALL Trap_Points(0.0d0,4.0d0*DASIN(1.0d0),nptsl,xleg,wleg)
!CALL W_Exp(xleg,nptsl,ExpMat)   

Hmatp = CMPLX(0.0d0,0.0d0)
Nmatp = CMPLX(0.0d0,0.0d0)
Hprime = CMPLX(0.0d0,0.0d0)

PHmatp = CMPLX(0.0d0,0.0d0)
PNmatp = CMPLX(0.0d0,0.0d0)

SDHmat = CMPLX(0.0d0,0.0d0)
SDNmat = CMPLX(0.0d0,0.0d0)

SDPHmat = CMPLX(0.0d0,0.0d0)
SDPNmat = CMPLX(0.0d0,0.0d0)

hmp = js

pairMat = 0.0d0
DO i = 1, nsps
    pairMat(i,i) = parOP(i)
END DO

! ///////////////// CALCULATE N_ijk & H_ijk  /////////////////
sditi = 1
sditj = 1

!!$OMP parallel do private(alp,bet,gam,psdr,nsdr,ovlpp,ovlpn,rhopij,rhonij,vme,RotMat)
!
!  NOTE: this will be more efficient if we "unroll" the three alp, bet,gam loops into one long loop
!        for future work

DO alp = 1, int(2.d0*js)+1 !np = bigJmax
    DO bet = 1, int(js)+1
        DO gam = 1, int(2.d0*js)+1

            CALL Psi_New(alpha_i(alp),beta_j(bet),gamma_k(gam),RotMat)

            !Modified for strictly Gauss-Legendre quadrature
            !Using generated R matrix, Rotate the neutron and proton slater determinants  
            psdr = MATMUL(RotMat,psd(sditj,:,:))
            nsdr = MATMUL(RotMat,nsd(sditj,:,:))

            ! Generate neutron and proton overlaps (ovlp*) as well as the density matrices (rho*ij)
            CALL makerhoij(1,numprot,psdr,psd(sditi,:,:),ovlpp,rhopij)
            CALL makerhoij(2,numneut,nsdr,nsd(sditi,:,:),ovlpn,rhonij)

            N_ijk(alp,bet,gam) = ovlpp*ovlpn

            allOvlpp(1,alp,bet,gam) = ovlpp
            allOvlpn(1,alp,bet,gam) = ovlpn

            psdr = psd(sditj,:,:)
            nsdr = nsd(sditj,:,:)
            IF (.NOT. isNorm) THEN
                CALL TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
                H_ijk(alp,bet,gam) = vme
            END IF

            IF (ParityTest) THEN
                psdr = MATMUL(RotMat,psdr)
                nsdr = MATMUL(RotMat,nsdr)

                psdr = MATMUL(pairMat,psdr)
                nsdr = MATMUL(pairMat,nsdr)

                CALL makerhoij(1,numprot,psdr,psd(sditi,:,:),ovlpp,rhopij)
                CALL makerhoij(2,numneut,nsdr,nsd(sditi,:,:),ovlpn,rhonij)

                PN_ijk(alp,bet,gam) = ovlpp*ovlpn

                allOvlpp(2,alp,bet,gam) = ovlpp
                allOvlpn(2,alp,bet,gam) = ovlpn
                IF (.NOT. isNorm) THEN
                    CALL TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
                    PH_ijk(alp,bet,gam) = vme
                END IF
           END IF
        END DO
    END DO
END DO
!!$OMP end parallel do

END SUBROUTINE Projection_with_Parity

!=========================================================
!
!  creates a vector with the parities of each s.p. state
!
! OUTPUT:  parOp(i) = +/-1 for each s.p. state i
!
! CALLED BY
!
SUBROUTINE MakeParityOp(parOP)

USE system_parameters
USE spstate
USE psis

IMPLICIT NONE

REAL (KIND = 8), DIMENSION(nsps), INTENT(OUT) :: parOP
INTEGER :: i

parOP = 0.0d0
DO i = 1, nsps
    parOP(i) = (-1.0d0)**(spsqn(i)%l)
END DO

END SUBROUTINE MakeParityOp


!=========================================================
!
! original generalized eigensolution routines
! using diagonalization
!
!  INPUT:
!    ndim = dimension of matrices
!    nparity = # of parity states considered
!   Nmat(i,j) = ndim x ndim norm matrix
!   Hmat(i,j)= ndim x ndim hamiltonian matrix
!   tol = real variable, tolerance for eliminating zeros
!
! OUTPUT
!   posdef(ipar) = if positive definite states found (for parity ipar)
!   nf = # of states found
!   pevals(i=1,nf) = energy eigenvalues found
!   normvals(i=1,ndim)
!
! subroutines called:
!
!  zheev  -- lapack routine for diagonalizing a hermitian, complex*16 matrix
!
subroutine geneigsolverdiag(ndim,nparity,ipar,NMat,HMat,tol,posdef,nf,pevals,normvals)
	use errortests
	implicit none
	
	interface
		subroutine zhegvdiag(jobz,np,ndim,A,B,tol,aa,bb,vec,nf,info)  ! INTERFACE
			implicit none
		!------------ INPUTS -----------------	
		    CHARACTER          JOBZ
		    INTEGER            NP,NDIM
			complex*16    ::  a(np,ndim),b(np,ndim)
			real :: tol
		!------------ OUTPUTS ------------------
			real(kind=8)   :: aa(ndim),bb(ndim)
			complex*16     :: vec(np,ndim)
			integer :: nf
			integer :: info
		end subroutine zhegvdiag ! INTERFACE
	end interface	
	integer (kind=4) :: ndim, nparity   ! dimension of matrices

    integer ipar
	COMPLEX (KIND = 8):: Hmat(ndim,ndim),NMat(ndim,ndim)   !Input hamiltonian and norm matrices
	
	LOGICAL, INTENT(INOUT) :: posdef(2)
	REAL :: tol
	
	INTEGER (KIND = 8), INTENT(OUT) :: nf
	REAL (KIND = 8), INTENT(OUT) :: pevals(nparity,ndim)
	real (kind=8) :: normvals(nparity,ndim)

	integer (kind=4) :: nftmp	
		
	
	LOGICAL :: invFlag = .true.
	integer :: i,it,k,jtt
	integer :: info,lwork,pairP,ii,jj
	integer loop1

	real (kind = 8), ALLOCATABLE :: e(:),normies(:)
	complex*16, allocatable :: vec(:,:)

    !write(*,*) nf,np
	allocate(e(ndim),normies(ndim),vec(ndim,ndim))
	CALL zhegvdiag('n',ndim,ndim,Hmat,Nmat,tol,e,normies,vec,nftmp,info)
	if(info/=0)then
		print*,' Problem running zhegvdiag '
		return
	end if
	nf = nftmp
	pevals(ipar,:)=e(:)
	normvals(ipar,:)=normies(:)
	return
	
end	subroutine geneigsolverdiag
!=========================================================
!
!  zhegvdiag
!
!  routine for solving generalized eigenvalue problem A vec = aa B vec
!  by diagonalizing B
!
!  called lapack subroutine ZHEEV
!
! INPUT: jobz = 'N' or 'n' for eigenvalues only, 'V' or 'v' for eigenvectors
!        uplo = 'U' or 'L' for stored in upper or lower triangular
!        np   = declared dimension of the arrays
!        ndim = dimension used by the arrays
!        A(:,:),B(:,:) Hermitian, double complex input matrices
!        tol : real valued tolerance for discarding eigenvalues of B
! OUTPUT:
!      aa(:) : final eigenvalues from A vec = aa B vec
!      bb(:) : eigenvalues of matrx B; used as check
!      vec(:,:) : if jobz = 'v', filled afterwards with eigenvectors
!      nf : # of final eigenvalues kept ( = # of eigenvalues of B > tol )
!      info: if =0, then results okay
!

subroutine zhegvdiag(jobz,np,ndim,A,B,tol,aa,bb,vec,nf,info)

	implicit none
	interface
	   SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,  INFO ) ! INTERFACE
	        CHARACTER          JOBZ, UPLO
	        INTEGER            INFO, LDA, LWORK, N
	        DOUBLE PRECISION   RWORK( * ), W( * )
	        COMPLEX*16         A( LDA, * ), WORK( * )
		end subroutine zheev      ! INTERFACE
    end interface	
!------------ INPUTS -----------------	
    CHARACTER          JOBZ
    INTEGER            NP,NDIM
	complex*16    ::  a(np,ndim),b(np,ndim)
	real :: tol
!------------ OUTPUTS ------------------
	
	real(kind=8)   :: aa(ndim),bb(ndim)
	complex*16     :: vec(np,ndim)
	integer :: nf
	integer :: info

!-------- INTERMEDIATE----------------
    complex*16, allocatable :: c(:,:)
    real (kind = 8),ALLOCATABLE :: rwork(:)
    COMPLEX (kind = 8), ALLOCATABLE :: work(:)
	integer :: lwork
	integer :: i,j,k,l,ii,jj
	complex*16 :: ztmp
	real(kind=8) :: facti,factj
	
!.......... FIRST STEP..... DIAGONALIZE B..................

    allocate(c(ndim,ndim))
	do i = 1,ndim
		do j = 1,ndim
			c(i,j)=b(i,j)
		end do
	end do
    lwork = 2*ndim - 1
    ALLOCATE(work(lwork))
    ALLOCATE(rwork(3*ndim - 2))
!	print*,' norm '
!	do i = 1,ndim
!		write(6,'(9f10.4)')(real(c(i,j)),j=1,ndim)
!	end do
	CALL zheev('v','u',ndim,c,ndim,bb,work,lwork,rwork,info)
	if(info/=0)then
		deallocate(c,work,rwork)
		return
	end if
	nf = 0
	do i = 1,ndim
		if(bb(i) > tol)then
			nf = nf + 1
		end if
	end do
!	print*,nf,' states '
	if(nf ==0)then
		deallocate(c,work,rwork)
		return
	end if
!................. NOW CONSTRUCT NEW MATRIX..... V = 
    ii = 0
    do i = 1,ndim
		if(bb(i)> tol)then
			ii = ii +1
			facti = 1.d0/dsqrt(bb(i))
		else
			cycle
		end if
		jj = 0
		do j = 1 ,ndim
			if(bb(j)> tol)then
				jj = jj +1
				factj = 1.d0/dsqrt(bb(j))
			else
				cycle
			end if
			ztmp = (0.d0,0.d0)
			do k = 1,ndim
				do l = 1,ndim
  				    ztmp = ztmp + dconjg(c(k,i))*a(k,l)*c(l,j)
				end do
			end do
			vec(ii,jj)= ztmp*facti*factj
			
		end do		
	end do
!	print*,' H '
!	do i = 1,ndim
!		write(6,'(9f10.4)')(real(a(i,j)),j=1,ndim)
!	end do
!	print*,' transformed '
!	do i = 1,nf
!		write(6,'(9f10.4)')(real(vec(i,j)),j=1,nf)
!	end do
!... NOW DIAGONALIZE....
	CALL zheev('v','u',nf,vec,np,aa,work,lwork,rwork,info)
!    print*,nf,' eigenvalues ',(aa(i),i=1,nf)
!.... TRANSFORM EIGENVECTORS IF NEEDED...TO BE ADDED LATER...	
	
	deallocate(c,work,rwork)
	
	return
	
end subroutine zhegvdiag
!=========================================================
!
! original generalized eigensolution routines
! using Cholesky
!
!  INPUT:
!    ndim = dimension of matrices
!    nparity = # of parity states considered
!   Nmat(i,j) = ndim x ndim norm matrix
!   Hmat(i,j)= ndim x ndim hamiltonian matrix
!   tol = real variable, tolerance for eliminating zeros
!
! OUTPUT
!   posdef(ipar) = if positive definite states found (for parity ipar)
!   nf = # of states found
!   pevals(i=1,nf) = eigenvalues found
!
! subroutines called:
!  choleskymaster -- finds cholesky decomposition
!  choleskyinvert -- inverts cholesky decomposition and applys to a matrix
!  zheev  -- lapack routine for diagonalizing a hermitian, complex*16 matrix
!
subroutine geneigsolverCholesky(ndim,nparity,ipar,NMat,HMat,tol,posdef,nf,pevals)
	use errortests
	implicit none
	
	interface
    
		  subroutine choleskyinvert(np,n,nf,A,B,Linv,skip)   ! INTERFACE
	      implicit none
	      integer NP,N 	! dimension of all array
	      integer nf        ! final dimension after skipping
	      complex(kind = 8) :: A(NP,NP),B(Np,NP),Linv(NP,np)
	      integer ::           skip(NP)
	      end subroutine choleskyinvert  ! INTERFACE
	
	      subroutine choleskymaster(np,n,A,tol,invflag,L, Linv, skip,posdef)  ! INTERFACE
	      implicit none
	      integer NP,N 	! dimension of all array
	      complex(kind = 8) :: A(NP,NP), L(NP,NP),Linv(NP,np)
	      integer ::           skip(NP)
	      real    ::           tol
	      logical ::           invflag, posdef
	      end subroutine choleskymaster ! INTERFACE
	  
		  SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,  INFO ) ! INTERFACE
		        CHARACTER          JOBZ, UPLO
		        INTEGER            INFO, LDA, LWORK, N
		        DOUBLE PRECISION   RWORK( * ), W( * )
		        COMPLEX*16         A( LDA, * ), WORK( * )
			end subroutine zheev      ! INTERFACE
	end interface	
	integer (kind=4) :: ndim, nparity   ! dimension of matrices

    integer ipar
	COMPLEX (KIND = 8):: Hmat(ndim,ndim),NMat(ndim,ndim)   !Input hamiltonian and norm matrices
	
	LOGICAL, INTENT(INOUT) :: posdef(2)
	REAL :: tol
	
	INTEGER (KIND = 8), INTENT(OUT) :: nf
	REAL (KIND = 8), INTENT(OUT) :: pevals(nparity,ndim)

	integer (kind=4) :: nftmp	
	
	COMPLEX (KIND = 8),allocatable :: Hprime(:,:)
	
	COMPLEX (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: L, Linv
	INTEGER, ALLOCATABLE, DIMENSION(:) :: skip
	
	LOGICAL :: invFlag = .true.
	integer :: i,it,k,jtt
	integer :: info,lwork,pairP,ii
	integer loop1

	real (kind = 8), ALLOCATABLE :: e(:)
	real (kind = 8),ALLOCATABLE :: rwork(:)
	COMPLEX (kind = 8), ALLOCATABLE :: work(:)

    ALLOCATE(L(ndim,ndim),Linv(ndim,ndim),skip(ndim))
	allocate(Hprime(ndim,ndim))

    CALL choleskymaster(ndim,ndim,Nmat,tol,invFlag,L,Linv,skip,posdef(ipar))
    IF (.NOT.diagp) THEN
        pevals(ipar,:) = 1002.0d0
        DEALLOCATE(L,Linv,skip)
        return
    END IF
    IF (posdef(ipar) .EQV. .true.) THEN
        CALL choleskyinvert(ndim,ndim,nftmp,Hmat,Hprime,Linv,skip) ! Returns H N^-1
		nf=int(nftmp,kind=8)
!        IF (nf > nkeep) THEN
!            nkeep = nf
!        END IF
    ELSE
        !WRITE(*,*) 'Matrix not positive-definite, skipping',posdef(pairOUT),pairOUT
        pevals(ipar,:) = 1000.0d0
        DEALLOCATE(L,Linv,skip,hprime)
        return
    ENDIF
    IF (nf == 0) THEN
        !WRITE(*,1099)
        DEALLOCATE(L,Linv,skip,hprime)
        return
    ENDIF
		
    lwork = 2*nf - 1
    ALLOCATE(work(lwork))
    ALLOCATE(rwork(3*nf - 2))
    ALLOCATE(e(nf))
    e = 0.0d0
    !write(*,*) nf,np
	CALL zheev('n','u',int(nf,kind=4),Hprime,ndim,e,work,lwork,rwork,info)
    DO loop1 = 1, nf
        !WRITE(*,*) (DBLE(Hprime(loop1,loop2)),loop2 = 1, nf)
        pevals(ipar,loop1) = e(loop1)
		
        !write(*,*) e(loop1)
    END DO
    DEALLOCATE(L,Linv,skip,e,work,rwork,Hprime)
	return
	
end	subroutine geneigsolverCholesky
!!=========================================================
