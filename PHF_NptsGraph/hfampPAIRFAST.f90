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
!    Js - User defined angular momentum (coupled - integer (or half-integer) value)
!    np - dimension of generated H and N matrixes (np = 2*floatJ + 1)
!    mm, mp - indices for generated matrix
!    hmm, hmp - actual values of m', m for wignerD matrix generation
!    uw - Generated WignerD matrix from user Js
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


SUBROUTINE Projection_with_Parity(jmin,jmax,isOdd,psd,nsd,tol,npts,nf,posdef,ParTest,ntrace, &
													problist,jlist,probsum,nftotal,jall,pallPair)

USE phf_vals
USE system_parameters
USE spstate
USE psis
USE errortests
USE jaggedArrayType

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

real (kind=8), allocatable, dimension(:,:) :: normvals
integer (kind=4) :: ilist
real (kind=8), allocatable, dimension(:,:) :: pevals
integer (kind=8), intent(out) :: nftotal
real (kind=8), intent(in) :: jmin, jmax
real (kind=8), dimension(int((jmax-jmin+1.0d0)*(jmin+jmax+1.0d0))), intent(inout) :: jall
real (kind=8), dimension(2,int((jmax-jmin+1.0d0)*(jmin+jmax+1.0d0))), intent(inout) :: pallPair
real(kind=8), intent(inout) :: problist(2,int(jmax-jmin)+1), jlist(int(jmax-jmin)+1), probsum
logical, intent(in) :: isOdd
integer (kind=8) :: intJ
real (kind=8) :: floatJ
INTEGER (KIND=8) :: np
COMPLEX (KIND = 8), INTENT(IN) :: psd(numsd,nsps,numprot), nsd(numsd,nsps,numneut)
LOGICAL, INTENT(IN) :: ParTest
LOGICAL, INTENT(INOUT) :: posdef(2)
INTEGER (KIND = 8), INTENT(OUT) :: nf
real (kind=8) :: ntrace(2)

integer (kind=4) :: nftmp
INTEGER (KIND = 8) :: sditi,sditj
INTEGER (KIND = 8) :: mp,mm,npts,nptsl,alp,bet,gam,errorflag,loop1,loop2,tt,hcheck,nkeep,itt
INTEGER, ALLOCATABLE, DIMENSION(:) :: skip
REAL :: tol
REAL (KIND = 8) :: a,b,hmp,hm,lam,pi,dx
REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: xleg,wleg,xlegb,wlegb
!REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:,:) :: xjac,wjac
COMPLEX (KIND = 8) :: uw,jerkap,jerkapn,jerkbp,jerkbpn,jerkgp,jerkgpn,ic
COMPLEX (KIND = 8) :: Pjerkap,Pjerkapn,Pjerkbp,Pjerkbpn,Pjerkgp,Pjerkgpn

type (jaggedArray), dimension(int(jmin):int(jmax)) :: phfA, phfB, phfG, NphfA, NphfB, NphfG
type (jaggedArray), dimension(int(jmin):int(jmax)) :: PphfA, PphfB, PphfG, PNphfA, PNphfB, PNphfG

COMPLEX (KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: RotMat, L, Linv,RotTesta,RotTestb,RotTestg
LOGICAL :: invFlag = .true.
PARAMETER(pi = 2.0d0*DASIN(1.0d0))
integer :: i,it,k,jtt

COMPLEX (KIND = 8),ALLOCATABLE :: psdr(:,:),nsdr(:,:),psdpr(:,:),nsdpr(:,:),ExpMat(:,:,:)
COMPLEX (KIND = 8),ALLOCATABLE :: rhopij(:,:),rhonij(:,:)

type (jaggedArray), dimension(int(jmin):int(jmax)) :: IntWig
type (jaggedArray), dimension(int(jmin):int(jmax)) :: Hmat, Nmat, Hprime
type (jaggedArray), dimension(int(jmin):int(jmax)) :: Hmatp, Hmatn, Nmatp, Nmatn
type (jaggedArray), dimension(int(jmin):int(jmax)) :: PHmatp, PHmatn, PNmatp, PNmatn

!Multiple Slater Determinant matrices
type (jaggedArray), dimension(int(jmin):int(jmax)) :: SDHmat, SDNmat, SDPHmat, SDPNmat
type (jaggedArray), dimension(int(jmin):int(jmax)) :: TSDHmat,TSDNmat
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


IF (ParTest) CALL MakeParityOp(parOP)

!npts = 35  ! set from outside
!tol = 0.1e0  ! set from outside
ALLOCATE(xleg(npts),wleg(npts))
ALLOCATE(xlegb(npts),wlegb(npts))
!ALLOCATE(xjac(np,np,npts),wjac(np,np,npts))

ALLOCATE(RotTesta(nsps,nsps),RotTestb(nsps,nsps),RotTestg(nsps,nsps))
ALLOCATE(RotMat(nsps,nsps),ExpMat(nsps,nsps,npts ))
ALLOCATE (rhopij(nsps,nsps),rhonij(nsps,nsps))

do intJ = int(jmin), int(jmax)
	np = 2*intJ+1
	if (isOdd) np = np + 1
	allocate(phfA(intJ)%MK(np,np),phfB(intJ)%MK(np,np),phfG(intJ)%MK(np,np))
	allocate(NphfA(intJ)%MK(np,np),NphfB(intJ)%MK(np,np),NphfG(intJ)%MK(np,np))
	allocate(PphfA(intJ)%MK(np,np),PphfB(intJ)%MK(np,np),PphfG(intJ)%MK(np,np))
	allocate(PNphfA(intJ)%MK(np,np),PNphfB(intJ)%MK(np,np),PNphfG(intJ)%MK(np,np))

	allocate(Hmatp(intJ)%MK(np,np),Hmatn(intJ)%MK(np,np),Nmatp(intJ)%MK(np,np),Nmatn(intJ)%MK(np,np))
	allocate(PHmatp(intJ)%MK(np,np),PHmatn(intJ)%MK(np,np),PNmatp(intJ)%MK(np,np),PNmatn(intJ)%MK(np,np))
	allocate(Hmat(intJ)%MK(np,np),Nmat(intJ)%MK(np,np),Hprime(intJ)%MK(numsd*np,numsd*np))

	allocate(IntWig(intJ)%MK(np,np))

	allocate(SDHmat(intJ)%MK(numsd*np,numsd*np),SDNmat(intJ)%MK(numsd*np,numsd*np))
	allocate(SDPHmat(intJ)%MK(numsd*np,numsd*np),SDPNmat(intJ)%MK(numsd*np,numsd*np))
	allocate(TSDHmat(intJ)%MK(numsd*np,numsd*np),TSDNmat(intJ)%MK(numsd*np,numsd*np))

	Hmatp(intJ)%MK = cmplx(0.0d0,0.0d0)
	Nmatp(intJ)%MK = cmplx(0.0d0,0.0d0)
	Hprime(intJ)%MK = cmplx(0.0d0,0.0d0)

	PHmatp(intJ)%MK = cmplx(0.0d0,0.0d0)
	PNmatp(intJ)%MK = cmplx(0.0d0,0.0d0)

	SDHmat(intJ)%MK = cmplx(0.0d0,0.0d0)
	SDNmat(intJ)%MK = cmplx(0.0d0,0.0d0)

	SDPHmat(intJ)%MK = cmplx(0.0d0,0.0d0)
	SDPNmat(intJ)%MK = cmplx(0.0d0,0.0d0)

end do

!Compute weights and abscissas for the alpha and gamma
!integrals using Gauss-Legendre quadrature

CALL gauleg(0.0d0,4.0d0*DASIN(1.0d0),xleg,wleg,npts)
CALL gauleg(0.0d0,2.0d0*DASIN(1.0d0),xlegb,wlegb,npts)
!CALL Trap_Points(0.0d0,4.0d0*DASIN(1.0d0),nptsl,xleg,wleg)
!CALL W_Exp(xleg,nptsl,ExpMat)   


! Main quadrature loop
pairMat = 0.0d0
DO i = 1, nsps
	pairMat(i,i) = parOP(i)
END DO
DO sditi = 1, numsd
	DO sditj = 1, numsd
		do intJ = int(jmin), int(jmax)
		    phfA(intJ)%MK =  CMPLX(0.0d0,0.0d0)
		    NphfA(intJ)%MK = CMPLX(0.0d0,0.0d0)
		    IF (ParTest) PphfA(intJ)%MK = CMPLX(0.0d0,0.0d0)
		    IF (ParTest) PNphfA(intJ)%MK = CMPLX(0.0d0,0.0d0)
		    !RotTesta = CMPLX(0.0d0,0.0d0)
		end do
        DO alp = 1, npts
			do intJ = int(jmin), int(jmax)
		        phfB(intJ)%MK = CMPLX(0.0d0,0.0d0)
		        NphfB(intJ)%MK = CMPLX(0.0d0,0.0d0)
		        IF (ParTest) PphfB(intJ)%MK = CMPLX(0.0d0,0.0d0)
		        IF (ParTest) PNphfB(intJ)%MK = CMPLX(0.0d0,0.0d0)
		        !RotTestb = CMPLX(0.0d0,0.0d0) 2.000   2.000
			end do
            DO bet = 1, npts
				do intJ = int(jmin), int(jmax)
		            phfG(intJ)%MK = CMPLX(0.0d0,0.0d0)
		            NphfG(intJ)%MK = CMPLX(0.0d0,0.0d0)
		            IF (ParTest) PphfG(intJ)%MK = CMPLX(0.0d0,0.0d0)
		            IF (ParTest) PNphfG(intJ)%MK = CMPLX(0.0d0,0.0d0)
		            !RotTestg = CMPLX(0.0d0,0.0d0)
				end do
                DO gam = 1, npts
                    CALL Psi_New(xleg(alp),xlegb(bet),xleg(gam),RotMat)

                    !	Modified for strictly Gauss-Legendre quadrature
                    !Using generated R matrix, Rotate the neutron and proton slater determinants
                    psdr = MATMUL(RotMat,psd(sditj,:,:))
                    nsdr = MATMUL(RotMat,nsd(sditj,:,:))

                    ! Generate neutron and proton overlaps (ovlp*) as well as the density matrices (rho*ij)
                    CALL makerhoij(1,numprot,psdr,psd(sditi,:,:),ovlpp,rhopij)
                    CALL makerhoij(2,numneut,nsdr,nsd(sditi,:,:),ovlpn,rhonij)

                    !Generate <H> = vme for each step in the integration
                    CALL TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)

					do intJ = int(jmin), int(jmax)
						floatJ = dble(intJ)
						if (isOdd) floatJ = floatJ + 0.5d0

						np = int(2.0d0*floatJ) + 1
						projfactor = (2*floatJ+1.d0)/8.d0/(pi)**2

		                call Wigner_d2(floatJ,int(np,kind=4),xleg(alp),xlegb(bet),xleg(gam),IntWig(intJ)%MK)

	                    ! build the integral for each m' and m in your user defined Wigner matrix
		                phfG(intJ)%MK = phfG(intJ)%MK + wleg(gam)*IntWig(intJ)%MK*vme*DSIN(xlegb(bet))*projfactor
		                NphfG(intJ)%MK = NphfG(intJ)%MK + wleg(gam)*IntWig(intJ)%MK*ovlpp*ovlpn*DSIN(xlegb(bet))*projfactor
					end do

                    !Recompute with Pairity operator on slater determinant
                    !    BEGIN PAIRITY CALCULATIONS
                    !Apply pairity operator
                    
                    psdr = psd(sditj,:,:)
                    nsdr = nsd(sditj,:,:)
                    
                    IF (ParTest) then
                    
                        psdr = MATMUL(RotMat,psdr)
                        nsdr = MATMUL(RotMat,nsdr)
                        
			          	psdr = MATMUL(pairMat,psdr)
			          	nsdr = MATMUL(pairMat,nsdr)
			          	
			          	! Hit it with the parity operator
			          	
			          	!FORALL(pairP = 1:nsps)
                        !    psdr(pairP,:) = pairOP(pairP)*psdr(pairP,:)
                        !    nsdr(pairP,:) = pairOP(pairP)*nsdr(pairP,:)
                        !END FORALL

              	        ! Generate neutron and proton overlaps (ovlp*) as well as the density matrices (roh*ij)
                        
                        !psdr = MATMUL(RotMat,psd(sditj,:,:))
                        !nsdr = MATMUL(RotMat,nsd(sditj,:,:))
                        
                        ! THEN you rotate
                        
                        !psdr = MATMUL(RotMat,psdr)
                        !nsdr = MATMUL(RotMat,nsdr)

                    	CALL makerhoij(1,numprot,psdr,psd(sditi,:,:),ovlpp,rhopij)
                    	CALL makerhoij(2,numneut,nsdr,nsd(sditi,:,:),ovlpn,rhonij)

                    	!Generate <H> = vme for each step in the integration
					!	print*,' parity transformed ',real(ovlpp*ovlpn)

                    	CALL TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)

						do intJ = int(jmin), int(jmax)		! MOVE TO MAIN LOOP, I.E. "vme2"
							floatJ = dble(intJ)
							if (isOdd) floatJ = floatJ + 0.5d0

							np = int(2.0d0*floatJ) + 1
							projfactor = (2*floatJ+1.d0)/8.d0/(pi)**2

				            call Wigner_d2(floatJ,int(np,kind=4),xleg(alp),xlegb(bet),xleg(gam),IntWig(intJ)%MK)

			                ! build the integral for each m' and m in your user defined Wigner matrix
		                	PphfG(intJ)%MK = PphfG(intJ)%MK + wleg(gam)*IntWig(intJ)%MK*vme*DSIN(xlegb(bet))*projfactor
		                	PNphfG(intJ)%MK = PNphfG(intJ)%MK + wleg(gam)*IntWig(intJ)%MK*ovlpp*ovlpn*DSIN(xlegb(bet))*projfactor
						end do
 	
		           END IF

                END DO
				do intJ = int(jmin), int(jmax)
		            phfB(intJ)%MK = phfB(intJ)%MK + wlegb(bet)*phfG(intJ)%MK
		            NphfB(intJ)%MK = NphfB(intJ)%MK + wlegb(bet)*NphfG(intJ)%MK
		            IF (ParTest) PphfB(intJ)%MK = PphfB(intJ)%MK + wlegb(bet)*PphfG(intJ)%MK
		            IF (ParTest) PNphfB(intJ)%MK = PNphfB(intJ)%MK + wlegb(bet)*PNphfG(intJ)%MK
		            !RotTestb = RotTestb + wlegb(bet)*RotTestg
				end do
            END DO
			do intJ = int(jmin), int(jmax)
		        phfA(intJ)%MK = phfA(intJ)%MK + wleg(alp)*phfB(intJ)%MK
		        NphfA(intJ)%MK = NphfA(intJ)%MK + wleg(alp)*NphfB(intJ)%MK
		        IF (ParTest) PphfA(intJ)%MK = PphfA(intJ)%MK + wleg(alp)*PphfB(intJ)%MK
		        IF (ParTest) PNphfA(intJ)%MK = PNphfA(intJ)%MK + wleg(alp)*PNphfB(intJ)%MK
	            !RotTesta = RotTesta + wleg(alp)*RotTestb
			end do
        END DO   ! alp
	do intJ = int(jmin), int(jmax)
		floatJ = dble(intJ)
		if (isOdd) floatJ = floatJ + 0.5d0
		np = int(2.0d0*floatJ) + 1

		Hmatp(intJ)%MK = phfA(intJ)%MK
		!phf(jloc)%h = Hmatp
		!Store H-matrix in appropriate location in larger Mult-SD matrix
		SDHmat(intJ)%MK((sditi-1)*np+1:sditi*np,(sditj-1)*np+1:sditj*np) = Hmatp(intJ)%MK

		Nmatp(intJ)%MK = NphfA(intJ)%MK
		!phf(jloc)%n = Nmatp
		!Store N-matrix in appropriate location in larger Mult-SD matrix
		SDNmat(intJ)%MK((sditi-1)*np+1:sditi*np,(sditj-1)*np+1:sditj*np) = Nmatp(intJ)%MK

		!Repeat with mixed-parity matrices
		IF (ParTest) THEN
			PHmatp(intJ)%MK = PphfA(intJ)%MK
			!phfp(jloc)%h = PHmatp
		  	SDPHmat(intJ)%MK((sditi-1)*np+1:sditi*np,(sditj-1)*np+1:sditj*np) = PHmatp(intJ)%MK

			PNmatp(intJ)%MK = PNphfA(intJ)%MK
			!phfp(jloc)%n = PNmatp
			SDPNmat(intJ)%MK((sditi-1)*np+1:sditi*np,(sditj-1)*np+1:sditj*np) = PNmatp(intJ)%MK
		END IF

	end do

	!End the multiple slater determinant loops
	END DO
END DO

!	TEST MATRICES
!
!testing = .false.
!if(testing)then
!write(*,*) "H Matrix"
!do ii = 1,np
!	write(*,*) SDHmat(ii,:)
!end do
!write(*,*)
!write(*,*) "PH Matrix"
!do ii = 1,np
!	write(*,*) SDPHmat(ii,:)
!end do
!write(*,*) 
!write(*,*) "N Matrix"
!do ii = 1,np
!	write(*,*) SDNmat(ii,:)
!end do
!write(*,*) 
!write(*,*) "PN Matrix"
!do ii = 1,np
!	write(*,*) SDPNmat(ii,:)
!end do
!write(*,*)
!write(*,*) "Energies ::"
!end if
!
!	Compute Trace of the Norm matrix (testing)
!
!trnorm = 0.0d0
!do ii = 1, np
!	trnorm = trnorm + REAL(Nmatp(ii,ii))
!end do
!
!write(*,*) "HF Contribution: ",trnorm

1099 format('No states projected for J =  ',1F8.1)
1001 FORMAT(I6,4X,F12.5,4X,F8.5,4X,F3.1,4x)
1002 FORMAT(I6,4x,2(F13.5),4x,F3.1)

nftotal = 0
probsum = 0.0d0
ilist = 0
allocate(normvals(npair,numsd*np))
do intJ = int(jmin), int(jmax)
	floatJ = dble(intJ)
	if (isOdd) floatJ = floatJ + 0.5d0
	np = int(2.0d0*floatJ) + 1
	nkeep = 0
	allocate(pevals(npair,numsd*np))
	pevals = 0.0d0
	ntrace = 0.0d0

	!Use Cholesky decomposition to produce H N^-1
	DO  pairOUT = 1, npair
		TSDHmat(intJ)%MK = CMPLX(0.0d0,0.0d0)
		TSDNmat(intJ)%MK = CMPLX(0.0d0,0.0d0)
		IF (ParTest) THEN
			TSDHmat(intJ)%MK = SDHmat(intJ)%MK + (-1.0d0)**(pairOUT+1)*SDPHmat(intJ)%MK
			TSDNmat(intJ)%MK = SDNmat(intJ)%MK + (-1.0d0)**(pairOUT+1)*SDPNmat(intJ)%MK
		ELSE
			TSDHmat(intJ)%MK = SDHmat(intJ)%MK
			TSDNmat(intJ)%MK = SDNmat(intJ)%MK
		END IF
		!nf = numsd*np

		!Computing trace of norm matrix		
		DO itt = 1,np
		    ntrace(pairout) = ntrace(pairout) + real(TsdNmat(intJ)%MK(itt,itt))
		END DO    

	!............. ROUTINE FOR SOLVING GENERALIZED PROBLEM...........
		call geneigsolverdiag(int(numsd*np,kind=4),npair,pairOUT,TSDNMat(intJ)%MK,TSDHMat(intJ)%MK,tol,posdef,nf,pevals,normvals)
		    IF (nf > nkeep) THEN
		        nkeep = nf
		    END IF
	END DO
	nf = nkeep

	! Post-Eigensolver -> write to screen and record for later printing
	ilist=ilist+1
	if(ParTest)then
	    problist(:,ilist)=ntrace(:)/float(numsd)
		probsum = probsum+ntrace(1)+ntrace(2)
	else
	    problist(1,ilist)=ntrace(1)/float(numsd)
		probsum = probsum + ntrace(1)/float(numsd)
	end if
	jlist(ilist)=floatJ

	IF ((posdef(1) .EQV. .TRUE.).OR.(posdef(2) .EQV. .TRUE.)) THEN
		DO ii = 1, nf
			nftotal = nftotal + 1
			IF (ParTest) THEN
				WRITE(*,1002) nftotal,pevals(1,ii),pevals(2,ii),floatJ
			ELSE
				WRITE(*,1001) nftotal,pevals(1,ii),floatJ
			END IF
			jall(nftotal) = floatJ
			pallPair(1,nftotal) = pevals(1,ii)
			IF (ParTest) pallPair(2,nftotal) = pevals(2,ii)
		END DO
	END IF

    if(nf<1)then
		write(*,*)' No states for J = ',floatJ
	end if
	deallocate(pevals)
end do

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
