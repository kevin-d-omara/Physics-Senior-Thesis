! ///////////////// VERSION CONTROL /////////////////
!	
! Latest Version 1.1.8.4
!
! note: this document was written with a tab width of 4 spaces
!
!	Changelog
!   1.0.0   10/02/15    : finalized paper outline of easyProjection, began implementation
!   1.0.1   10/03/15    : implemented initialization and N_ijk overlap matrix population
!   1.0.2   10/04/15    : implemented outer M and K loops and "populateZeta" subroutine for Z_Mi & Z_Kk
!                           - vetted both loop indices and spot checked populateZeta routine
!                           - success of populateZeta implies success of alpha_i & gamma_k
!                           - note: I have a hand written page with more detailed testing notes
!                         & added blueprints with testing notes
!                         & added populateN_jMK subroutine
!   1.0.3   10/04/15    : added inner J' & J indices, Jmin, sizeOfNorm, and floatJ indices
!                           - vetted both inner loop indice & floatJ indices
!                         & added 2 subroutines:
!                           - putElementIn_N_Jprime_MK
!                           - putElementIn_deltaJpJ
!   1.0.4   10/04/15    : added subroutine ZGESV -> inverts and solves for the norm (N_J_MK)
!                           - briefly checked "populateN_ijk": each element abs(value) <= 1
!   1.0.5   10/04/15    : added method to add trace components as they are calculated
!                           - IF FINALLY WORKS!!!!!!!!
!                           - now to see if these results match those from the integration method
!   1.0.6   10/08/15    : added relevant Hamiltonian matrices (mirror of Norm matrices)
!   1.0.7   10/09/15    : added variable definitions to the preamble (dictionary)
!                         re-worded subroutines involving the Norm and Hamiltonian
!                           -populateNH_jMK & putElementIn_NH_Jprime_MK
!   1.0.8   10/10/15    : added jaggedArray type in module jaggedArrayType (in phf_modules_full.f)
!                         changed norm & Hamiltonian solutions (N_J_MK & H_J_MK) to jaggedArray types
!                         added subroutine extractPiecemealElements
!                         added subroutine mirrorSolution
!   1.0.9   10/10/15    : added geneigsolverdiag to find eigenvalues of Hamiltonian
!                           -correct number of states found, but the energies are off
!   1.0.10  10/12/15    : the norm and Hamiltonian are complex hermitian matrices, therefore the
!                         symmetric values must be complex conjugates
!                           - modified subroutine mirrorSolution to include "dconjg" for complex conjugates
!   1.0.11  10/16/15    : cleaned up code, moved eigensolving to subroutine "EigenSolverPackage"
!                           - updated dictionary
!   1.1.0   10/16/15    : found the error! ZGESV "destroys" the deltaJpJ term, which was being used
!                         subsequently to find the Hamiltonian
!                           - added tempDeltaJpJ to rectify this issue
!   1.1.1   10/22/15    : began work on including parity -> updated EigenSolverPackage
!   1.1.2   10/22/15    : finished including parity
!   1.1.3   11/03/15    : added subroutine for user input of tolerance and numberOfPoints (npts)
!   1.1.4   11/03/15    : fully implemented change to a user entered number of points (npts)
!   1.1.5   11/12/15    : added feature to track how long the program takes (elapsedTime)
!   1.1.5.1 11/12/15    : implemented "write results to screen" (J_WriteResults)
!                           - need to tidy up and add definitions
!   1.1.6   11/15/15    : added sum of trace(H)
!                           - combined with trace(N) and moved to single subroutine "addToTrace"
!   1.1.6.1 11/15/15    : output to screen and optional file now mirrors phf -> less confusion this way
!   1.1.6.2 11/15/15    : fixed issue with normSum and hamSum parity
!   1.1.6.3 11/18/15    : fixed issue with sum of trace of norm/Hamiltonian (wasn't looping large enough)
!   1.1.6.5 12/03/15    : updated "J_Write_Results" to inclue the frac(J) for the Hamiltonian
!   1.1.7   12/03/15    : fixed issue with improper parity trace values
!                           - was adding parity term after trace calculations: <|> +- <Psi|RP|Psi>
!                           - removed printTrace and addToTrace subroutines
!	1.1.7.1	03/11/16	: quick patch -> sum of Hamiltonian no longer reports sum of Norm
!   1.1.8.0 04/11/16    : begin Independent Norm => initial unentanglement of Norm & Hamiltonian
!   1.1.8.1 04/12/16    : added extractNorm & printNorm
!   1.1.8.2 04/12/16    : independent norm COMPLETE! only missing "autoChoose" feature
!   1.1.8.3 04/13/16    : changed terminology: ephf.x --> lamp.x (linear algebra angular momentum projection)
!                           - replaced npts with bigJmax
!   1.1.8.4 04/13/16    : updated phf_applyh.f and phf_tbme.f
!                           - from phf_v10 to phf_v13
!                           - properly ordered loops (i.e. column-major; row fastest changing index for speed)
!
! Todo:
!   - add title info at first run which includes version number
!   - add "autoChooser"
!   - replace tabs with spaces
!
! Note: '!!' means not in dictionary yet (also not in subroutine description)
!
!  POTENTIAL PROBLEMS:
!    - wigner_d function -> M, K possibly flipped
!    - subroutine extractPiecemealElements & mirrorSolution not heavily tested
!    - all work for extracting eigenvalues is spread thin (i.e. little testing, little understanding b.c. I'm borrowing from J. Staker)
!    - replace Tabs with Spaces; note: some already replaced, but not all
!
! /////////////////// BLUEPRINTS ////////////////////
!
! initialize
! enter Jmax (& determine bigJmax & numOfJ)
!
! allocate(N_ijk(bigJ,bigJ,numJ))
! allocate(alp(bigJ),bet(numJ),gam(bigJ))
! allocate(Z_Mi(bigJ),Z_Kk(bigJ))
! allocate(N_jMK(numJ))
! allocate(N_J_MK(sizeOfNorm))
!
! build alp(i), bet(j), gam(k)
! build N_ijk & H_ijk
!
! do M = -Jmax, +Jmax                               !indices vetted
!   build Z_Mi                                      !vetted(spot): populateZeta(Z_Mi)
!   do K = -Jmax, M !M>K                            !indices vetted
!     build Z_Kk
!     build N_jMK & H_jMK                           !briefly check for abs(value) <= 1
!
!     Jmin = max(abs(floatM),abs(floatK))
!     sizeOfNorm = int(Jmax - Jmin) + 1
!     allocate(N^J'_MK(sizeOfNorm))
!     allocate(deltaJpJ(sizeOfNorm,sizeOfNorm))
!
!     do J' = 1, sizeOfNorm                         !indices vetted
!       build N^J'_MK & H_^J'_MK !1 element each pass
!       do J = 1, sizeOfNorm                        !indices vetted
!         build deltaJpJ         !1 row each pass
!       end J'
!     end J
!
!     INVERT (eq. 22)
!     deallocate(N^J'_MK,deltaJpJ)
!   end K
! end M
!
! ///////////////////////////////////////////////////
!****************************************************
!
! 'Easy' projection of angular momentum
!
!    Dr. Calvin Johnson
!    Kevin O'Mara
!
!****************************************************
!
! Main inversion routine: (dictionary)
!
! //initialize//
!   psd, nsd = proton/neutron slater determinants
!   psdf, nsdf = ^^ with a third dimension to allow multiple slater determinants (i.e. numsd > 1)
!   ParityTest = logical, .true. if even/odd parity
!
! //inputJmax//
!   Jmax = maximum J value for projection [user input]
!   bigJmax = 2*Jmax + 1
!   isOdd = logical flag for odd A nuclides
!   numOfJ = # of beta values -> int(Jmax + 1.)
!
! //tolerance & number of points//
!   tolerance = cuttoff point for eliminating zeroes in the Norm matrix (during zhegvdiag)
!
! //allocate (pre-body)//
!   N_ijk, H_ijk = norm/Hamiltonian overlap matrix [eq. 13]
!   allOvlpp(2,i,j,k) = ovlpp => dimension 1: 1=even 2=odd
!   allOvlpn(2,i,j,k) = ovlpn => dimension 1: 1=even 2=odd
!   allRhopij(2,i,j,k,:,:) = rhopij => dimension 1: 1=even 2=odd
!   allRhonij(2,i,j,k,:,:) = rhopij => dimension 1: 1=even 2=odd
!   alpha_i, gamma_k = vector containing angles / overlap points [before eq. 17]
!   beta_j = vector containing angles / overlap points [eq. 20]
!   Z_Mi, Z_Kk = inverse alpha/gamma solution [eq. 17]
!   N_jMK, H_jMK = first intermediate norm/hamiltonian vector for a fixed value of M & K [eq. 18]
!   N_J_MK, H_J_MK = norm/Hamiltonian jagged array containing full solution with indices (J,M,K) [eq. 22]
!   nlist = length of problist, hamlist, and jlist; int(Jmax-Jmin)+1
!
! //build N_ijk & H_ijk//
!   posdefP = if positive definite states found [logical]
!   dbleInt = function call to convert an Integer(kind=4) to Integer(king=8)
!
! //main body//
!   intM, intK = matrix index for M & K
!   floatMK = *shifted* real value for M or K (shifted to accomodate for negative M or K values)
!            ie: intM = -2, -1,  0, +1, +2
!              floatM = +1, +2, +3, +4, +5
!   Jmin = minimum J value -> max(abs(floatM),abs(floatK))
!   sizeOfNorm = size of the Norm and Hamiltonian matrices -> int(Jmax - Jmin) + 1
!   intJprime, intJ = vector index for Jprime & J
!   floatJprime, floatJ = *shifted* real value for corresponding intJprime/intJ
!   N_Jprime_MK, H_Jprime_MK = second intermediate norm/Hamiltonian vector for a fixed M & K [eq. 21]
!   deltaJpJ = singular value decomposition (SVD) object for a fixed M & K [eq. 19]
!   IPIVOT = [not used]
!   INFO = contains error flag i.e. if INFO /= 0 then there is an issue
!
! //trace//
!   normTrace = vector containing the trace of each N^J_MK matrix (fixed J, indices M & K)
!   hamTrace = vector containing the trace of each H^J_MK matrix (fixed J, indices M & K)
!
subroutine J_Main

	use phf_vals
	use system_parameters
	use spstate
	use sporbit
	use psis
	use errortests
	use jaggedArrayType

	implicit none

	interface
		subroutine allocateSlaterDet(psd,nsd,psdf,nsdf,ParityTest)
			use system_parameters
			use spstate
			implicit none
			complex (kind = 8),allocatable, intent(out) :: psd(:,:),nsd(:,:)
			complex (kind = 8),allocatable, intent(out) :: psdf(:,:,:),nsdf(:,:,:)
			logical, intent(out) :: ParityTest
		end subroutine allocateSlaterDet

		subroutine inputJmax(Jmax,bigJmax,isOdd,numOfJ)
			use system_parameters
			implicit none
			real (kind=4), intent(out) :: Jmax
			integer, intent(out) :: bigJmax
			logical, intent(out) :: isOdd
			integer, intent(out) :: numOfJ
		end subroutine inputJmax

		subroutine inputTolerance(tolerance)
			implicit none
			real (kind=4), intent(out) :: tolerance
		end subroutine inputTolerance

		subroutine populateAlpBetGam(bigJmax,numOfJ,alpha_i,beta_j,gamma_k)
			implicit none
			integer, intent(in) :: bigJmax, numOfJ
			real (kind=8), allocatable, intent(out) :: alpha_i(:), gamma_k(:), beta_j(:)
		end subroutine populateAlpBetGam

		subroutine Projection_with_Parity(js,np,psd,nsd,tol,ParTest,alpha_i,beta_j, &
            gamma_k,N_ijk,H_ijk,PN_ijk,PH_ijk,allOvlpp,allOvlpn,allRhopij,allRhonij,isNorm)
			use phf_vals
			use system_parameters
			use spstate
			use psis
			use errortests
			implicit none
			integer (kind = 8), intent(in) :: np
			real (kind = 8), intent(in) :: js
			complex (kind = 8), intent(in) :: psd(numsd,nsps,numprot), nsd(numsd,nsps,numneut)
			real (kind=4) :: tol
			logical, intent(in) :: ParTest
			real (kind=8), allocatable, intent(in) :: alpha_i(:), gamma_k(:), beta_j(:)
			complex (kind = 8), allocatable, intent(inout)	:: N_ijk(:,:,:), H_ijk(:,:,:), PN_ijk(:,:,:), PH_ijk(:,:,:)
            complex (kind = 8), allocatable, intent(inout) :: allOvlpp(:,:,:,:), allOvlpn(:,:,:,:)
            complex (kind = 8), allocatable, intent(inout) :: allRhopij(:,:,:,:,:,:), allRhonij(:,:,:,:,:,:)
            logical, intent(in) :: isNorm
	   end subroutine Projection_with_Parity

		subroutine populateZeta(intMK,floatMK,bigJmax,angle_ik,zeta)
			implicit none
			integer (kind=4), intent(in) :: intMK, bigJmax
			real (kind=4), intent(in) :: floatMK
			real (kind=8), allocatable, intent(in) :: angle_ik(:)
			complex (kind=8), allocatable, intent(inout) :: zeta(:)
		end subroutine populateZeta

		subroutine populateNH_jMK(Jmax,bigJmax,numOfJ,Z_Mi,Z_Kk,NH_ijk,NH_jMK)
			implicit none
			real (kind=4), intent(in) :: Jmax
			integer (kind=4), intent(in) :: bigJmax, numOfJ
			complex(kind=8), allocatable, intent(in) :: Z_Mi(:), Z_Kk(:)
			complex (kind = 8), allocatable, intent(in)	:: NH_ijk(:,:,:)
			complex (kind = 8), allocatable, intent(inout) :: NH_jMK(:)
		end subroutine populateNH_jMK

		subroutine putElementIn_NH_Jprime_MK(intJprime,floatJprime,numOfJ,floatM,floatK,beta_j,NH_jMK,NH_Jprime_MK)
			implicit none
			integer (kind=4), intent(in) :: intJprime, numOfJ
			real (kind=4), intent(in) :: floatJprime, floatM, floatK
			real (kind=8), allocatable, intent(in) :: beta_j(:)	!kind=8 -> kind=4 on call because wignerD
			complex (kind = 8), allocatable, intent(in) :: NH_jMK(:)
			complex (kind = 8), allocatable, intent(inout) :: NH_Jprime_MK(:)
		end subroutine putElementIn_NH_Jprime_MK

		subroutine putElementIn_deltaJpJ(intJprime,intJ,floatJprime,floatJ,floatM,floatK,beta_j,numOfJ,deltaJpJ)
			implicit none
			integer (kind=4), intent(in) :: intJprime, intJ, numOfJ
			real (kind=4), intent(in) :: floatJprime, floatJ, floatM, floatK
			real (kind=8), allocatable, intent(in) :: beta_j(:)
			complex (kind = 8), allocatable, intent(inout) :: deltaJpJ(:,:)
		end subroutine putElementIn_deltaJpJ

		subroutine errorZGESV(INFO)
			implicit none
			integer, intent(in) :: INFO
		end subroutine errorZGESV

		subroutine extractPiecemealElements(numOfJ,sizeOfNorm,intM,intK,invertedSolution,NH_J_MK)
			use jaggedArrayType
			implicit none
			integer (kind=4), intent(in) :: numOfJ, sizeOfNorm, intM, intK
			complex (kind = 8), allocatable, intent(in) :: invertedSolution(:)
			type (jaggedArray), allocatable, intent(inout) :: NH_J_MK(:)
		end subroutine extractPiecemealElements

		subroutine mirrorSolution(NH_J_MK)
			use jaggedArrayType
			implicit none
			type (jaggedArray), allocatable, intent(inout) :: NH_J_MK(:)
		end subroutine mirrorSolution

        subroutine extractNorm(numOfJ,isOdd,N_J_MK,PN_J_MK,ParityTest,nlist,problist)
        	use phf_vals
	        use jaggedArrayType
	        implicit none
	        integer (kind=4), intent(in) :: numOfJ
	        logical, intent(in) :: isOdd
	        type (jaggedArray), allocatable, intent(in) :: N_J_MK(:), PN_J_MK(:)
	        logical, intent(in) :: ParityTest
	        integer (kind=4), intent(in) :: nlist
	        real (kind=8), intent(out) :: problist(2,nlist)
        end subroutine extractNorm

        subroutine printNorm(numOfJ,isOdd,ParityTest,nlist,problist,printTo)
	        implicit none
	        integer (kind=4), intent(in) :: numOfJ
	        logical, intent(in) :: isOdd
	        logical, intent(in) :: ParityTest
	        integer (kind=4), intent(in) :: nlist
	        real (kind=8), intent(out) :: problist(2,nlist)
            integer (kind=4) :: printTo
        end subroutine printNorm

		subroutine EigenSolverPackage(numOfJ,isOdd,N_J_MK,H_J_MK,PN_J_MK,PH_J_MK,posdefP,ParityTest,tolerance, &
										npmax,nftotal,jall,pallPair,nlist,jlist,normSum,hamSum,problist,hamlist)
			use phf_vals
			use jaggedArrayType
			implicit none
			integer (kind=4), intent(in) :: numOfJ
			logical, intent(in) :: isOdd
			type (jaggedArray), allocatable, intent(in) :: N_J_MK(:), H_J_MK(:), PN_J_MK(:), PH_J_MK(:)
			logical, intent(inout) :: posdefP(2)
			logical, intent(in) :: ParityTest
			real (kind=4), intent(in) :: tolerance
			integer (kind=8), intent(in) :: npmax!!
			integer (kind=4), intent(in) :: nlist!!
			integer (kind=8), intent(out) :: nftotal!!
			real (kind=8), intent(inout) :: jall(npmax), pallPair(2,npmax)!!
			real (kind=8), intent(inout) :: jlist(nlist)!!
			complex (kind = 8), intent(out) :: normSum, hamSum
			real (kind=8), intent(out) :: problist(2,nlist), hamlist(2,nlist)!!
		end subroutine EigenSolverPackage
	end interface

    !initialize
	complex (kind=8), allocatable :: psd(:,:),nsd(:,:)
	complex (kind=8), allocatable :: psdf(:,:,:),nsdf(:,:,:)
	logical :: ParityTest
	logical :: test1 = .true.

    !inputJmax
	real (kind=4) :: Jmax
	integer (kind=4) :: bigJmax !2*Jmax+1
	logical :: isOdd
	integer (kind=4) :: numOfJ	!int(Jmax+1.)

	!tolerance & number of points
	integer (kind=4) :: clockStart, clockEnd, clockRate
	real (kind=4) :: tolerance, elapsedTime

    !allocate (pre-body)
	complex (kind=8), allocatable :: N_ijk(:,:,:), H_ijk(:,:,:), PN_ijk(:,:,:), PH_ijk(:,:,:)
    complex (kind=8), allocatable :: allOvlpp(:,:,:,:), allOvlpn(:,:,:,:)
    complex (kind = 8), allocatable :: allRhopij(:,:,:,:,:,:), allRhonij(:,:,:,:,:,:)
	real (kind=8), allocatable :: alpha_i(:), gamma_k(:), beta_j(:)
	complex(kind=8), allocatable :: Z_Mi(:), Z_Kk(:)
	complex (kind=8), allocatable :: N_jMK(:), H_jMK(:), PN_jMK(:), PH_jMK(:)
	type (jaggedArray), allocatable :: N_J_MK(:), H_J_MK(:), PN_J_MK(:), PH_J_MK(:)
	integer (kind=4) :: tempBigJ
	integer (kind=8) :: npmax!!
	integer (kind=4) :: nlist!!
	real (kind=8), allocatable :: jall(:), pallPair(:,:)!!
	real (kind=8), allocatable :: jlist(:), problist(:,:), hamlist(:,:)!!
	integer (kind=8) :: nftotal!!doesn't go here

	! build N_ijk & H_ijk
	logical :: posdefP(2) = .true.
	integer (kind=8) :: dbleInt !function

    ! main body
	integer (kind=4) :: intM, intK
	real (kind=4) :: floatM, floatK
	real (kind=4) :: Jmin
	integer (kind=4) :: sizeOfNorm
	integer (kind=4) :: intJprime, intJ
	real (kind=4) :: floatJprime, floatJ
	complex (kind=8), allocatable :: N_Jprime_MK(:), H_Jprime_MK(:), PN_Jprime_MK(:), PH_Jprime_MK(:)
	complex (kind=8), allocatable :: deltaJpJ(:,:), tempDeltaJpJ(:,:), PNtempDeltaJpJ(:,:), PHtempDeltaJpJ(:,:)
	integer (kind=4), allocatable :: IPIVOT(:)
	integer (kind=4) :: INFO

    ! prompt user
    character (len=1) :: choice
    logical :: autoChoose, finalNorm
    integer :: iostatus
    character (len=100) :: filename

	! trace
	complex (kind=8), allocatable :: normTrace(:), PnormTrace(:), hamTrace(:), PhamTrace(:)
	complex (kind = 8) :: normSum, hamSum

! ##################################################################################################

	!initialize
	call testset
	call allocateSlaterDet(psd,nsd,psdf,nsdf,ParityTest)
	call PairLog(ParityTest)
	call inputTolerance(tolerance)

    autoChoose = .false.
    finalNorm = .false.
    normLoop: do
	    !input Jmax
        if (autoChoose) then
!            call autoChooser()
            print *, "autoChooser (under construction)"
        else
    	    call inputJmax(Jmax,bigJmax,isOdd,numOfJ)
        end if

	    call system_clock(COUNT_RATE = clockRate)
	    call system_clock(COUNT = clockStart)
	    !allocate (pre-body)
	    allocate(N_ijk(bigJmax,numOfJ,bigJmax),H_ijk(bigJmax,numOfJ,bigJmax))
	    allocate(PN_ijk(bigJmax,numOfJ,bigJmax),PH_ijk(bigJmax,numOfJ,bigJmax))
        allocate(allOvlpp(2,bigJmax,numOfJ,bigJmax),allOvlpn(2,bigJmax,numOfJ,bigJmax))
        allocate(allRhopij(2,bigJmax,numOfJ,bigJmax,nsps,nsps),allRhonij(2,bigJmax,numOfJ,bigJmax,nsps,nsps))
	    allocate(alpha_i(bigJmax),beta_j(numOfJ),gamma_k(bigJmax))
	    allocate(Z_Mi(bigJmax),Z_Kk(bigJmax))
	    allocate(N_jMK(numOfJ), H_jMK(numOfJ))
	    if (ParityTest) allocate(PN_jMK(numOfJ), PH_jMK(numOfJ))
	    allocate(normTrace(numOfJ),PnormTrace(numOfJ),hamTrace(numOfJ),PhamTrace(numOfJ))
	    allocate(N_J_MK(numOfJ),H_J_MK(numOFJ))
	    if (ParityTest) allocate(PN_J_MK(numOfJ),PH_J_MK(numOFJ))
	    do intJ =1, numOfJ !allocate "jagged" sizes to norm and Hamiltonian
    !		tempBigJ = int(2.*intJ)-1
		    tempBigJ = 2*intJ-1
		    if (isOdd) tempBigJ = tempBigJ + 1
		    allocate(N_J_MK(intJ)%MK(tempBigJ,tempBigJ))
		    allocate(H_J_MK(intJ)%MK(tempBigJ,tempBigJ))
		    if (ParityTest) then
			    allocate(PN_J_MK(intJ)%MK(tempBigJ,tempBigJ))
			    allocate(PH_J_MK(intJ)%MK(tempBigJ,tempBigJ))
		    end if
	    end do
        Jmin = 0.0
        if (isOdd) Jmin = 0.5
	    npmax = int((Jmax-Jmin+1.0d0)*(Jmin+Jmax+1.0d0))!!
	    nlist = int(Jmax-Jmin)+1
	    allocate(jall(npmax),pallPair(2,npmax),jlist(nlist),problist(2,nlist),hamlist(2,nlist))!!
	
	    call populateAlpBetGam(bigJmax,numOfJ,alpha_i,beta_j,gamma_k)

	    ! build N_ijk & H_ijk
	    if (.not.ParityTest) posdefP(2) = .false. !necessary(?)
	    call Projection_with_Parity(dble(Jmax),dbleInt(bigJmax),psdf,nsdf,tolerance, & !populates N_ijk, H_ijk
                ParityTest,alpha_i,beta_j,gamma_k,N_ijk,H_ijk,PN_ijk,PH_ijk,allOvlpp,allOvlpn,allRhopij,allRhonij,.true.)

	    ! main body (Norm)
	    do intM = 1, bigJmax !-Jmax to Jmax
		    floatM = float(intM) - (Jmax + 1.)
		    call populateZeta(intM,floatM,bigJmax,alpha_i,Z_Mi) !Z_Mi [eq. 17]

		    do intK = 1, intM !-Jmax to M; M>=J
			    floatK = float(intK) - (Jmax + 1.)
			    call populateZeta(intK,floatK,bigJmax,gamma_k,Z_Kk) !Z_Kk [eq. 17]
			    call populateNH_jMK(Jmax,bigJmax,numOfJ,Z_Mi,Z_Kk,N_ijk,N_jMK) !first intermediate Norm [eq. 18]
			    if (ParityTest) then
				    call populateNH_jMK(Jmax,bigJmax,numOfJ,Z_Mi,Z_Kk,PN_ijk,PN_jMK)
			    end if

			    Jmin = max(abs(floatM),abs(floatK))
			    sizeOfNorm = int(Jmax - Jmin) + 1
			    allocate(N_Jprime_MK(sizeOfNorm))
			    if (ParityTest) allocate(PN_Jprime_MK(sizeOfNorm))
			    allocate(deltaJpJ(sizeOfNorm,sizeOfNorm),tempDeltaJpJ(sizeOfNorm,sizeOfNorm))
			    if (ParityTest) allocate(PNtempDeltaJpJ(sizeOfNorm,sizeOfNorm))
			    allocate(IPIVOT(sizeOfNorm))

			    do intJprime = 1, sizeOfNorm ! Jmin to Jmax
				    floatJprime = float(intJprime) + Jmin - 1.															!second intermediate
				    call putElementIn_NH_Jprime_MK(intJprime,floatJprime,numOfJ,floatM,floatK,beta_j,N_jMK,N_Jprime_MK) !Norm [eq. 21]
				    if (ParityTest) then
					    call putElementIn_NH_Jprime_MK(intJprime,floatJprime,numOfJ,floatM,floatK,beta_j,PN_jMK,PN_Jprime_MK)
				    end if

				    do intJ = 1, sizeOfNorm ! Jmin to Jmax
					    floatJ = float(intJ) + Jmin - 1.
					    call putElementIn_deltaJpJ(intJprime,intJ,floatJprime,floatJ,floatM,floatK,beta_j,numOfJ,deltaJpJ) !deltaJpJ [eq. 19]
				    end do
			    end do

			    ! Invert to find the Norm and Hamiltonian matrices, N_J_MK & H_J_MK
    !			tempDeltaJpJ = deltaJpJ
			    if (ParityTest) then
				    PNtempDeltaJPJ = deltaJpJ
			    end if

			    ! Solve Norm [eq 22] (held in N_Jprime_MK) & build full solution N_J_MK
			    call ZGESV(sizeOfNorm,1,deltaJpJ,sizeOfNorm,IPIVOT,N_Jprime_MK,sizeOfNorm,INFO)
			    call errorZGESV(INFO)
			    call extractPiecemealElements(numOfJ,sizeOfNorm,intM,intK,N_Jprime_MK,N_J_MK)

			    if (ParityTest) then
				    call ZGESV(sizeOfNorm,1,PNtempDeltaJPJ,sizeOfNorm,IPIVOT,PN_Jprime_MK,sizeOfNorm,INFO)
				    call errorZGESV(INFO)
				    call extractPiecemealElements(numOfJ,sizeOfNorm,intM,intK,PN_Jprime_MK,PN_J_MK)
			    end if

			    deallocate(N_Jprime_MK,deltaJpJ,tempDeltaJpJ,IPIVOT)
			    if (ParityTest) deallocate(PN_Jprime_MK,PNtempDeltaJpJ)
		    end do
	    end do

	    call mirrorSolution(N_J_MK)
	    if (ParityTest) then
		    call mirrorSolution(PN_J_MK)
	    end if

        if (finalNorm) exit normLoop
        call extractNorm(numOfJ,isOdd,N_J_MK,PN_J_MK,ParityTest,nlist,problist)
        call printNorm(numOfJ,isOdd,ParityTest,nlist,problist,6)

        call system_clock(COUNT = clockEnd)
        elapsedTime = (real(clockEnd) - real(clockStart))/real(clockRate)
        print '(/,A,F8.2,A,/)', "Total norm calculation time: ", elapsedTime, " seconds"
        ! end main body (Norm)

        ! prompt user
        promptLoop: do
            print *, "Select a choice"
            print *, "(C) Continue with Hamiltonian calculation"
            print *, "(R) Re-calculate norm w/different maxJ"
            print *, "(P) Print norm results to a file"
            print *, "(X) Exit"
            read(*,*) choice
            normChoice: select case(choice)
                case('c','C')
                    exit normLoop
                case('r','R')
                    print *, "(A) Automatically determine lowest maxJ (under construction)"
                    print *, "(M) Manually chose"
                    read(*,*) choice
                    autoChoose = .false.
                    if ((choice == 'a').OR.(choice == 'A')) autoChoose = .true.
                    print *, "(F) Final norm calculation (i.e. no prompt, continue to Hamiltonian afterwards)"
                    print *, "(N) Not final calculation (i.e. display result and prompt again)"
                    read(*,*) choice
                    finalNorm = .false.
                    if ((choice == 'f').OR.(choice == 'F')) finalNorm = .true.

                    deallocate(N_ijk,H_ijk)
                    deallocate(PN_ijk,PH_ijk)
                    deallocate(allOvlpp,allOvlpn)
                    deallocate(allRhopij,allRhonij)
                    deallocate(alpha_i,beta_j,gamma_k)
                    deallocate(Z_Mi,Z_Kk)
                    deallocate(N_jMK, H_jMK)
                    if (ParityTest) deallocate(PN_jMK, PH_jMK)
                    deallocate(normTrace,PnormTrace,hamTrace,PhamTrace)
                    deallocate(N_J_MK,H_J_MK)
                    if (ParityTest) deallocate(PN_J_MK,PH_J_MK)
                    deallocate(jall,pallPair,jlist,problist,hamlist)!!

                    exit promptLoop

                case('p','P')
                    do
	                    write(*,*) 'Enter file name (without extention -- .dat added): '
	                    read(*,*) filename
	                    filename = TRIM(filename)//'.dat'
	                    open(unit = 66, file = filename, status = 'NEW',iostat = iostatus)
	                    if (iostatus > 0) then
		                    write(*,*) "File already exists: Overwrite (o), append (a), new file name (n)?"
		                    read(*,*) choice
		                    do
			                    if ((choice == 'o').or.(choice == 'O')) then
				                    close(unit = 66)
				                    open(unit = 66, file = filename, status = 'REPLACE')
				                    exit
			                    elseif ((choice == 'a').OR.(choice == 'A')) then
				                    close(UNIT = 66)
				                    open(UNIT = 66,file = filename, status = 'OLD', position = 'APPEND')
				                    exit
			                    elseif ((choice == 'n').OR.(choice == 'N')) then
				                    close(unit = 66)
				                    exit
			                    else
				                    write(*,*) 'Incorrect choice.  Please select again.'
			                    end if
		                    end do
		                    if ((choice.ne.'n').and.(choice.ne.'N')) exit
	                    else
		                    exit
	                    end if
                    end do

                    call printNorm(numOfJ,isOdd,ParityTest,nlist,problist,66)
                    print '(/,A,A,/)', 'Data written to:', filename
                case('x','X')
                    return
                case default
                    print *, "Invalid choice."
            end select normChoice
        end do promptLoop
    end do normLoop

	call system_clock(COUNT_RATE = clockRate)
	call system_clock(COUNT = clockStart)
	! main body (Hamiltonian)
	! build H_ijk (i.e. vme) (bottleneck)
	if (.not.ParityTest) posdefP(2) = .false. !necessary(?)
	call Projection_with_Parity(dble(Jmax),dbleInt(bigJmax),psdf,nsdf,tolerance, & !populates N_ijk, H_ijk
            ParityTest,alpha_i,beta_j,gamma_k,N_ijk,H_ijk,PN_ijk,PH_ijk,allOvlpp,allOvlpn,allRhopij,allRhonij,.false.)

	do intM = 1, bigJmax !-Jmax to Jmax
		floatM = float(intM) - (Jmax + 1.)
		call populateZeta(intM,floatM,bigJmax,alpha_i,Z_Mi) !Z_Mi [eq. 17]

		do intK = 1, intM !-Jmax to M; M>=J
			floatK = float(intK) - (Jmax + 1.)
			call populateZeta(intK,floatK,bigJmax,gamma_k,Z_Kk) !Z_Kk [eq. 17]
			call populateNH_jMK(Jmax,bigJmax,numOfJ,Z_Mi,Z_Kk,H_ijk,H_jMK) !first intermediate Hamiltonian [eq. 18]
			if (ParityTest) then
				call populateNH_jMK(Jmax,bigJmax,numOfJ,Z_Mi,Z_Kk,PH_ijk,PH_jMK)
			end if

			Jmin = max(abs(floatM),abs(floatK))
			sizeOfNorm = int(Jmax - Jmin) + 1
			allocate(H_Jprime_MK(sizeOfNorm))
			if (ParityTest) allocate(PH_Jprime_MK(sizeOfNorm))
			allocate(deltaJpJ(sizeOfNorm,sizeOfNorm),tempDeltaJpJ(sizeOfNorm,sizeOfNorm))
			if (ParityTest) allocate(PHtempDeltaJpJ(sizeOfNorm,sizeOfNorm))
			allocate(IPIVOT(sizeOfNorm))

			do intJprime = 1, sizeOfNorm ! Jmin to Jmax
				floatJprime = float(intJprime) + Jmin - 1.															!second intermediate
				call putElementIn_NH_Jprime_MK(intJprime,floatJprime,numOfJ,floatM,floatK,beta_j,H_jMK,H_Jprime_MK) !Hamiltonian [eq. 21]
				if (ParityTest) then
					call putElementIn_NH_Jprime_MK(intJprime,floatJprime,numOfJ,floatM,floatK,beta_j,PH_jMK,PH_Jprime_MK)
				end if

				do intJ = 1, sizeOfNorm ! Jmin to Jmax
					floatJ = float(intJ) + Jmin - 1.
					call putElementIn_deltaJpJ(intJprime,intJ,floatJprime,floatJ,floatM,floatK,beta_j,numOfJ,deltaJpJ) !deltaJpJ [eq. 19]
				end do
			end do

			! Invert to find the Norm and Hamiltonian matrices, N_J_MK & H_J_MK
!			tempDeltaJpJ = deltaJpJ
			if (ParityTest) then
				PHtempDeltaJpJ = deltaJpJ
			end if

			! Solve Hamiltonian [eq 22] (held in H_Jprime_MK) & build full solution H_J_MK
			call ZGESV(sizeOfNorm,1,deltaJpJ,sizeOfNorm,IPIVOT,H_Jprime_MK,sizeOfNorm,INFO)
			call errorZGESV(INFO)
			call extractPiecemealElements(numOfJ,sizeOfNorm,intM,intK,H_Jprime_MK,H_J_MK)

			if (ParityTest) then
				call ZGESV(sizeOfNorm,1,PHtempDeltaJPJ,sizeOfNorm,IPIVOT,PH_Jprime_MK,sizeOfNorm,INFO)
				call errorZGESV(INFO)
				call extractPiecemealElements(numOfJ,sizeOfNorm,intM,intK,PH_Jprime_MK,PH_J_MK)
			end if

			deallocate(H_Jprime_MK,deltaJpJ,tempDeltaJpJ,IPIVOT)
			if (ParityTest) deallocate(PH_Jprime_MK,PHtempDeltaJpJ)
		end do
	end do

	call mirrorSolution(H_J_MK)
	if (ParityTest) then
		call mirrorSolution(PH_J_MK)
	end if
    ! end main body (Hamiltonian)

	call EigenSolverPackage(numOfJ,isOdd,N_J_MK,H_J_MK,PN_J_MK,PH_J_MK,posdefP,ParityTest,tolerance, &
								npmax,nftotal,jall,pallPair,nlist,jlist,normSum,hamSum,problist,hamlist)

	call system_clock(COUNT = clockEnd)
	elapsedTime = (real(clockEnd) - real(clockStart))/real(clockRate)
	print '(/,A,F8.2,A)', "Total projection time: ", elapsedTime, " seconds"
	print *, ''
	print *,' Sum of norms = ', dble(normSum)
	print *,' Sum of trace(H) = ', dble(hamSum)

	call J_WriteResults(tolerance,npmax,nftotal,jall,pallPair,nlist,jlist,problist,hamlist,ParityTest)!!
! npmax?
! jall?
! pallPair? (+misspelled)
! nlist?
! jlist?
! problist?

return; end subroutine J_Main
!==============================================================================
!
! allocate the slater determinants psd / nsd & psdf / nsdf
!
! OUTPUT:
!	psd, nsd = proton/neutron slater determinants
!	psdf, nsdf = ^^ with a third dimension for multiple slater deterimants (i.e. numsd > 1)
!	ParityTest = checks for proper parity [?]
!
! SUBROUTINES CALLED:
!	GetSherpaSD = master routine for retrieving slater determinant written out by SHERPA
!	PairLog = checks for parity and flags PairLog appropriately
!
subroutine allocateSlaterDet(psd,nsd,psdf,nsdf,ParityTest)

	use system_parameters
	use spstate

	implicit none
!............OUTPUT......................
	complex (kind = 8),allocatable, intent(out) :: psd(:,:),nsd(:,:)
	complex (kind = 8),allocatable, intent(out) :: psdf(:,:,:),nsdf(:,:,:)
	logical, intent(out) :: ParityTest

!............INTERNAL....................
	integer :: ii

!----------------- ALLOCATE SLATER DETERMINANTS -------------------
	allocate (psd(nsps,numprot),nsd(nsps,numneut))
	numsd = 1				! why is this hard-coded?
	allocate (psdf(numsd,nsps,numprot),nsdf(numsd,nsps,numneut))
	do ii = 1,numsd
		call GetSherpaSD(psd,nsd)
		psdf(ii,:,:) = psd
		nsdf(ii,:,:) = nsd
	end do

return; end subroutine allocateSlaterDet
!==============================================================================
!
! user enters a value for Jmax -> the highest value to project
! also determines bigJmax, isOdd, and numOfJ
!
! OUTPUT:
!	Jmax = for projection [user input]
!	bigJmax = 2*Jmax + 1
!	isOdd = logical flag for odd A nuclides
!	numOfJ = # of beta values -> int(Jmax + 1.)
!
subroutine inputJmax(Jmax,bigJmax,isOdd,numOfJ)

	use system_parameters

	implicit none
!............OUTPUT......................
	real (kind=4), intent(out) :: Jmax
	integer, intent(out) :: bigJmax
	logical, intent(out) :: isOdd
	integer, intent(out) :: numOfJ

!............INTERNAL....................
	integer :: AA
	real (kind=8) :: check
	character (len = 12) :: jstring

!---------------------- USER ENTERS "Jmax" ------------------------
	AA = numprot + numneut
	if(mod(AA,2) == 0)then
		check = 0.0d0; jstring = 'integer'; isOdd = .false.
	else
		check = 1.0d0; jstring = 'half-integer'; isOdd = .true.
	endif

	write(*,*) 'Enter J-max to project (',trim(jstring),' value):'
	do
		read(*,*) Jmax
		if(mod(2.0d0*Jmax,2.0d0) == check)then
			exit
		else
			write(*,*) 'Incorrect j.  A = ', AA
			write(*,*) 'Expecting ', trim(jstring), ' values for J'
			write(*,*) 'Please re-enter J-max.'
		end if
	end do

!--------------------- SET "numOfJ" & "bigJmax" -----------------------
	numOfJ = int(Jmax + 1.)
	bigJmax = int(2.*Jmax)+1
!	if (test1) write(*,*) "bigJmax = ", bigJmax

return; end subroutine inputJmax
!==============================================================================
!
! user enters values for the Norm matrix tolerance and the number of points for inversion.
!
! OUTPUT:
!	tolerance = cuttoff point for eliminating zeroes in the Norm matrix (during zhegvdiag)
!
subroutine inputTolerance(tolerance)

	implicit none
!............OUTPUT......................
	real (kind=4), intent(out) :: tolerance

!---------------------- USER ENTERS VALUES ------------------------
	print *, ' Enter tolerance for norm (typical = 0.01) '
	read(*,*) tolerance

return; end subroutine inputTolerance
!==============================================================================
!
! [equation 17 (before)] and [equation 20]
! populates alpha_i, beta_j, and gamma_k vectors
!
!  INPUT:
!	bigJmax = 2*Jmax+1; number of alpha and gamma points in the norm/Hamiltonian overlap 
!
! OUTPUT:
!	alpha_i, gamma_k = vector containing angles / overlap points [before eq. 17]
!	beta_j = vector contanining angles / overlap points [eq. 20]
!
subroutine populateAlpBetGam(bigJmax,numOfJ,alpha_i,beta_j,gamma_k)

	implicit none
!........... INPUT........................	
	integer, intent(in) :: bigJmax, numOfJ

!............OUTPUT......................
	real (kind=8), allocatable, intent(out) :: alpha_i(:), beta_j(:), gamma_k(:)

!............INTERNAL....................
	real (kind=4), parameter :: pi = acos(-1.)
	integer :: ii

!-------------------- POPULATE ALPHA(i) & GAMMA(k) ----------------
	allocate(alpha_i(1:bigJmax)); allocate(gamma_k(1:bigJmax))

	do ii = 1, bigJmax
		alpha_i(ii) = (float(ii)-1.)*(2.*pi)/float(bigJmax)
		gamma_k(ii) = alpha_i(ii)
	!	if (test1) write(*,'(A,I2.1,A,F6.3)') 'alpha(',ii,') = ', alpha_i(ii)
	end do

!---------------------- POPULATE BETA(j) --------------------------
	allocate(beta_j(1:numOfJ))

	!if (test1) write(*,'(A,I2.1)') 'numOfJ = ', numOfJ
	do ii = 1, numOfJ
		beta_j(ii) = (float(ii) - .5)*pi/float(numOfJ)
	!	if (test1) write(*,'(A,I2.1,A,F7.4)') 'Beta ',ii,' = ', beta_j(ii)
	end do

return; end subroutine populateAlpBetGam
!==============================================================================
!
! takes INTEGER(kind=4) and returns INTEGER(kind=8)
!
function dbleInt(inputInt)
	implicit none
	integer (kind=4), intent(in) :: inputInt
	integer (kind=8) :: dbleInt

	dbleInt = inputInt

return; end function
!==============================================================================
!
! [equation 17]
! takes in a zeta vector (i.e. Z_Mi or Z_Kk) and populates it
!
!  INPUT:
!	intMK = vector index (M or K)
!	floatMK = *shifted* real value for M or K (shifted to accomodate for negative M or K values)
!	bigJmax = 2*Jmax+1; number of alpha and gamma points in the norm/Hamiltonian overlap 
!	angle_ik = alpha_i or gamma_k vector containing angles / overlap points [before eq. 17]
!
! OUTPUT:
!	zeta = Z_Mi or Z_Kk vector -> specific to the input M or K
!
subroutine populateZeta(intMK,floatMK,bigJmax,angle_ik,zeta)

	implicit none
!........... INPUT........................	
	integer (kind=4), intent(in) :: intMK, bigJmax
	real (kind=4), intent(in) :: floatMK
	real (kind=8), allocatable, intent(in) :: angle_ik(:)

!............OUTPUT......................
	complex (kind=8), allocatable, intent(inout) :: zeta(:)

!............INTERNAL....................
	integer (kind=4) :: ik

!------------------------- POPULATE ZETA --------------------------
	do ik = 1, bigJmax
		zeta(ik) = dcmplx(1.d0/dble(bigJmax),0.d0) * exp(-(0.d0,1.d0)*dble(floatMK)*angle_ik(ik))
!		if (test1) write(*,'(A,I2.1,A,I2.1,A,F7.4,2X,F7.4)') 'zeta(',intMK,',',ik,') = ', zeta(ik)
	end do

return; end subroutine populateZeta
!==============================================================================
!
! [equation 18]
! populates N_jMK or H_jMK vector for a fixed value of M & K
!
!  INPUT:
!	Jmax = for projection
!	bigJmax = 2*Jmax+1; number of alpha and gamma points in the norm/Hamiltonian overlap 
!	numOfJ = # of beta values
!	Z_Mi, Z_Kk = inverse alpha/gamma solution [eq. 17]
!	NH_ijk = norm or Hamiltonian overlap matrix (N_ijk or H_ijk) [eq. 13]
!
! OUTPUT:
!	NH_jMK = first intermediate norm/Hamiltonian vector for a fixed value of M & K
!
subroutine populateNH_jMK(Jmax,bigJmax,numOfJ,Z_Mi,Z_Kk,NH_ijk,NH_jMK)

	implicit none
!........... INPUT........................
	real (kind=4), intent(in) :: Jmax
	integer (kind=4), intent(in) :: bigJmax, numOfJ
	complex(kind=8), allocatable, intent(in) :: Z_Mi(:), Z_Kk(:)
	complex (kind = 8), allocatable, intent(in)	:: NH_ijk(:,:,:)

!............OUTPUT......................
	complex (kind = 8), allocatable, intent(inout) :: NH_jMK(:)

!............INTERNAL....................
	integer (kind=4) :: ii, jj, kk

!------------------------ POPULATE NH_jMK --------------------------
NH_jMK = (0.d0,0.d0)
do jj = 1, numOfJ
	do ii = 1, bigJmax
		do kk = 1, bigJmax
			NH_jMK(jj) = NH_jMK(jj) + Z_Mi(ii)*Z_Kk(kk)*NH_ijk(ii,jj,kk)
		end do
	end do
end do

return; end subroutine populateNH_jMK
!==============================================================================
!
! [equation 21]
! calculates current element for N_Jprime_MK or H_Jprime_MK with a fixed Jprime, M, and K
!
!  INPUT:
!	intJprime = index for Jprime
!	floatJprime = shifted real value for corresponding intJprime
!	numOfJ = # of beta values
!	floatM, floatK = shifted real value for M & K
!	beta_j = vector containing beta angles / overlap oints [eq. 20]
!	NH_jMK = first intermediate norm or Hamiltonian vector for a fixed value of M & K [eq. 18]
!
! OUTPUT:
!	NH_Jprime_MK(intJprime) = second intermediate norm/Hamiltonain vector -> single element (intJprime)
!
subroutine putElementIn_NH_Jprime_MK(intJprime,floatJprime,numOfJ,floatM,floatK,beta_j,NH_jMK,NH_Jprime_MK)

	implicit none

	interface
		function wigner_d(xjj,xmp,xm,theta)
			implicit none
			real(kind = 4) :: wigner_d
			real, intent(in) :: xjj,xmp,xm
			real,intent(in) :: theta
		end function wigner_d
	end interface

!........... INPUT........................
	integer (kind=4), intent(in) :: intJprime, numOfJ
	real (kind=4), intent(in) :: floatJprime, floatM, floatK
	real (kind=8), allocatable, intent(in) :: beta_j(:)
	complex (kind = 8), allocatable, intent(in) :: NH_jMK(:)

!............OUTPUT......................
	complex (kind = 8), allocatable, intent(inout) :: NH_Jprime_MK(:)

!............INTERNAL....................
	integer (kind=4) :: jj

!----------- CALCULATE SINGLE ELEMENT FOR N_Jprime_MK -------------
	NH_Jprime_MK(intJprime) = (0.d0,0.d0)
	do jj = 1, numOfJ
		NH_Jprime_MK(intJprime) = NH_Jprime_MK(intJprime) + &
		wigner_d(floatJprime,floatM,floatK,real(beta_j(jj)))*NH_jMK(jj)
	end do

return; end subroutine putElementIn_NH_Jprime_MK
!==============================================================================
!
! [equation 19]
! calculates current element for deltaJpJ -> i.e. fixed Jprime, J, M, and K
!
!  INPUT:
!	intJprime, intJ = index for Jprime, J
!	floatJprime, floatJ = shifted real value for corresponding intJprime, intJ
!	floatM, floatK = shifted real value for M & K
!	beta_j = vector containing beta angles / overlap points [eq. 20]
!	numOfJ = # of beta values
!
! OUTPUT:
!	deltaJpJ(intJprime,intJ) = deltaJpJ matrix -> single element (intJprime,intJ)
!
subroutine putElementIn_deltaJpJ(intJprime,intJ,floatJprime,floatJ,floatM,floatK,beta_j,numOfJ,deltaJpJ)

	implicit none

	interface
		function wigner_d(xjj,xmp,xm,theta)
			implicit none
			real(kind = 4) :: wigner_d
			real, intent(in) :: xjj,xmp,xm
			real,intent(in) :: theta
		end function wigner_d
	end interface

!........... INPUT........................
	integer (kind=4), intent(in) :: intJprime, intJ, numOfJ
	real (kind=4), intent(in) :: floatJprime, floatJ, floatM, floatK
	real (kind=8), allocatable, intent(in) :: beta_j(:)

!............OUTPUT......................
	complex (kind = 8), allocatable, intent(inout) :: deltaJpJ(:,:)

!............INTERNAL....................
	integer (kind=4) :: jj
	real (kind = 4) :: tempSum ! kind=4 because of wigner_d function

!------------- CALCULATE SINGLE ELEMENT FOR delatJpJ --------------
	tempSum = 0.
	do jj = 1, numOfJ
		tempSum = tempSum +	wigner_d(floatJprime,floatM,floatK,real(beta_j(jj)))* &
							wigner_d(floatJ,floatM,floatK,real(beta_j(jj)))
	end do
	deltaJpJ(intJprime,intJ) = dcmplx(tempSum)

return; end subroutine putElementIn_deltaJpJ
!==============================================================================
!
!	Checks for errors in ZGESV subroutine.
!
!  INPUT:
!	INFO = contains error flag i.e. if INFO /= 0 then there is an issue
!
subroutine errorZGESV(INFO)

	implicit none
!........... INPUT........................	
	integer, intent(in) :: INFO

!----------------------- CHECK FOR ERORS --------------------------
if (INFO > 0) then
	write(*,*) 'ERROR, factorization completed, but the factor U is '
	write(*,*) 'exactly singular, so the solution could not be computed.'
else if (INFO < 0) then
	write(*,'(A,I3.1)') 'ERROR, illegal value at element ', INFO
else
	continue
end if

return; end subroutine errorZGESV
!==============================================================================
!
! takes each element from the inverted solution and places it into the corresponding
! location in the Norm or Hamiltonian matrix (full solution)
!
!  INPUT:
!	numOfJ = # of beta values -> int(Jmax + 1.)
!	sizeOfNorm = size of the Norm and Hamiltonian matrices -> int(Jmax - Jmin) + 1
!	intM, intK = matrix index for M & K
!	invertedSolution = norm/Hamiltonian solution (N_Jprime_MK / H_Jprime_MK) [eq. 22]
!
! OUTPUT:
!	NH_J_MK = full solution (all J) for the norm or Hamiltonian [eq. 22]
!		**note: full solution is being constructed piecemeal, therefore it is not *complete* until the mirror
!
subroutine extractPiecemealElements(numOfJ,sizeOfNorm,intM,intK,invertedSolution,NH_J_MK)

	use jaggedArrayType

	implicit none
!........... INPUT........................
	integer (kind=4), intent(in) :: numOfJ, sizeOfNorm, intM, intK
	complex (kind = 8), allocatable, intent(in) :: invertedSolution(:)

!............OUTPUT......................
	type (jaggedArray), allocatable, intent(inout) :: NH_J_MK(:)

!............INTERNAL....................
	integer (kind=4) :: offSetJ, offSetMK, tempM, tempK, jj

!--------------------- EXTRACT ELEMENTS ---------------------------
	offSetJ = numOfJ - sizeOfNorm
	offSetMK = 1 - sizeOfNorm
	tempM = intM + offSetMK
	tempK = intK + offSetMK
!	write(*,*) 'numOfJ = ', numOfJ, ' sizeOfNorm = ', sizeOfNorm
!	write(*,*) 'offSetMK = ', offSetMK
!	write(*,*) 'intM = ', intM, ' intK = ', intK
!	write(*,*) 'tempM = ', tempM, ' tempK = ', tempK
!	write(*,*) ''
	do jj = 1, sizeOfNorm !Jmin to Jmax
		NH_J_MK(jj+offSetJ)%MK(tempM,tempK) = invertedSolution(jj)
		tempM = tempM + 1
		tempK = tempK + 1
	end do

return; end subroutine extractPiecemealElements
!==============================================================================
!
! takes the norm or Hamiltonian and mirrors the values from the lower triangle into the upper triangle
!
!  INPUT / OUTPUT:
!	NH_J_MK = full solution (all J) for the norm or Hamiltonian [eq. 22]
!
subroutine mirrorSolution(NH_J_MK)

	use jaggedArrayType

	implicit none
!..........INPUT / OUTPUT................
	type (jaggedArray), allocatable, intent(inout) :: NH_J_MK(:)

!............INTERNAL....................
	integer (kind=4) :: intJ, col, row

!---------------------- MIRROR VALUES -----------------------------
do intJ = 1, ubound(NH_J_MK,1)
	do col = 1, ubound(NH_J_MK(intJ)%MK,1)
		do row = 1+col, ubound(NH_J_MK(intJ)%MK,1)
			NH_J_MK(intJ)%MK(col,row) = dconjg(NH_J_MK(intJ)%MK(row,col))
		end do
	end do
end do

return; end subroutine mirrorSolution
!==============================================================================
!
! takes the solved Norm and Hamiltonian and runs them through "geneigsolverdiag" to eigensolve.
! after, prints out results to the screen.
!
!  INPUT:
!	numOfJ = # of beta values -> int(Jmax + 1.)
!   isOdd = logical flag for odd A nuclides
!	N_J_MK, H_J_MK = full norm/Hamiltonian solution, held in a jagged array [eq. 22]
!	posdefP = if positive definite states found [logical]
!   ParityTest = logical, .true. if even/odd parity
!   tolerance = cuttoff point for eliminating zeroes in the Norm matrix (during zhegvdiag)
!   nlist = length of problist, hamlist, and jlist; int(Jmax-Jmin)+1
!
! OUTPUT:
!   problist(2,:) = fraction of original HF state in the norm for each J
!
!  INTERNAL:
!	ParityEvals(i=1,numberFound) = energy eigenvalues found
!	normvals = not used
!	parityOUT = counts parity state (i.e. 1 or 2)
!	npair = number of parity states (i.e. 1 or 2) [via phf_vals module]
!	TSDHmat = temporary Slater determinant Hamiltonian matrix (holds values from H_J_MK)
!	TSDNmat = temporary Slater determinant Norm matrix (holds values from N_J_MK)
!	nf = number of states found
!
subroutine EigenSolverPackage(numOfJ,isOdd,N_J_MK,H_J_MK,PN_J_MK,PH_J_MK,posdefP,ParityTest,tolerance, &
									npmax,nftotal,jall,pallPair,nlist,jlist,normSum,hamSum,problist,hamlist)

	use phf_vals
	use jaggedArrayType

	implicit none

	interface
		subroutine geneigsolverdiag(ndim,npair,ipar,NMat,HMat,tol,posdef,nf,pevals,normvals)
			use errortests
			implicit none	
			integer (kind=4) :: ndim, npair   ! dimension of matrices
			integer ipar
			COMPLEX (KIND = 8):: Hmat(ndim,ndim),NMat(ndim,ndim)   !Input hamiltonian and norm matrices
			LOGICAL, INTENT(INOUT) :: posdef(2)
			REAL :: tol
			INTEGER (KIND = 8), INTENT(OUT) :: nf
			REAL (KIND = 8), INTENT(OUT) :: pevals(npair,ndim)
			real (kind=8) :: normvals(npair,ndim)
		end subroutine geneigsolverdiag
	end interface

!........... INPUT........................
	integer (kind=4), intent(in) :: numOfJ
	logical, intent(in) :: isOdd
	type (jaggedArray), allocatable, intent(in) :: N_J_MK(:), H_J_MK(:), PN_J_MK(:), PH_J_MK(:)
	logical, intent(inout) :: posdefP(2) !output not used
	logical, intent(in) :: ParityTest
	real (kind=4), intent(in) :: tolerance
	integer (kind=8), intent(in) :: npmax!!
	integer (kind=4), intent(in) :: nlist!!

!............OUTPUT......................
	integer (kind=8), intent(out) :: nftotal!!
	real (kind=8), intent(inout) :: jall(npmax), pallPair(2,npmax)!!
	real (kind=8), intent(inout) :: jlist(nlist)!!
	complex (kind = 8), intent(out) :: normSum, hamSum
	real (kind=8), intent(out) :: problist(2,nlist), hamlist(2,nlist)!!

!............INTERNAL....................
	integer (kind=4) :: intJ, kk
	real (kind=4) :: floatJ
	real (kind=8), allocatable :: ParityEvals(:,:), normvals(:,:)
	integer (kind=4) :: parityOUT
	complex(kind=8), allocatable, dimension(:,:) :: TSDHmat, TSDNmat
	integer (kind=8) :: nf, nkeep
	integer (kind=4) :: stateFound

	nftotal = 0!!
	jall = 0!!
	pallPair = 0.0d0!!
	jlist = 0.0d0!!
	normSum = 0.0d0
	hamSum = 0.0d0
	problist = 0.0d0
	hamlist = 0.0d0

	if (ParityTest) then
		write(*,*) ' State  E(+ parity)  E(-parity) J'		
	else
		write(*,*) ' State    E       J'
	end if
	write(*,*) ' ------------------------------------'

	do intJ = 1, numOfJ
		floatJ = float(intJ) - 1
		if (isOdd) floatJ = floatJ + 0.5
		allocate(ParityEvals(npair,ubound(N_J_MK(intJ)%MK,1)))
		allocate(normvals(npair,ubound(N_J_MK(intJ)%MK,1)))
		allocate(TSDHmat(ubound(H_J_MK(intJ)%MK,1),ubound(H_J_MK(intJ)%MK,2)))
		allocate(TSDNmat(ubound(N_J_MK(intJ)%MK,1),ubound(N_J_MK(intJ)%MK,2)))
		nkeep = 0
		ParityEvals = 0.0d0

		do parityOUT = 1, npair
			TSDHmat = cmplx(0.0d0,0.0d0)
			TSDNmat = cmplx(0.0d0,0.0d0)
			if (ParityTest) then
				TSDHmat = H_J_MK(intJ)%MK + (-1.0d0)**(parityOUT+1)*PH_J_MK(intJ)%MK
				TSDNmat = N_J_MK(intJ)%MK + (-1.0d0)**(parityOUT+1)*PN_J_MK(intJ)%MK
			else
				TSDHmat = H_J_MK(intJ)%MK
				TSDNmat = N_J_MK(intJ)%MK
			end if

			do kk = 1, ubound(TSDNmat,1)
				problist(parityOUT,intJ) = problist(parityOUT,intJ) + dble(TSDNmat(kk,kk))
				hamlist(parityOUT,intJ) = hamlist(parityOUT,intJ) + dble(TSDHmat(kk,kk))
			end do

			call geneigsolverdiag(ubound(TSDHmat,1),npair,parityOUT,TSDNmat,TSDHmat,tolerance,posdefP,nf,ParityEvals,normvals)
			if (nf > nkeep) nkeep = nf
		end do

		normsum = normsum + problist(1,intJ)
		hamsum = hamsum + hamlist(1,intJ)
		if (parityTest) then
			normsum = normsum + problist(2,intJ)
			hamsum = hamsum + hamlist(2,intJ)
		end if

		jlist(intJ) = floatJ!!
		do stateFound = 1, int(nkeep)
			nftotal = nftotal + 1
			if (ParityTest) then
				write(*,1001) nftotal,ParityEvals(1,stateFound),ParityEvals(2,stateFound),floatJ
			else
				write(*,1002) nftotal,ParityEvals(1,stateFound),floatJ
			end if
			jall(nftotal) = floatJ!!
			pallPair(1,nftotal) = ParityEvals(1,stateFound)!!
			if (ParityTest) pallPair(2,nftotal) = ParityEvals(2,stateFound)!!
		end do
        if (nkeep<1) write(*,*)' No states for J = ',floatJ

		deallocate(ParityEvals,normvals,TSDHmat,TSDNmat)
	end do

	1001 format(I3,2X,2F12.5,2X,F4.1)
	1002 format(I3,2X,F12.5,2X,F4.1)
return; end subroutine EigenSolverPackage

!==============================================================================
!
! populates problist from the sovled Norm.
!
!  INPUT:
!	numOfJ = # of beta values -> int(Jmax + 1.)
!   isOdd = logical flag for odd A nuclides
!	N_J_MK = full norm solution, held in a jagged array [eq. 22]
!   ParityTest = logical, .true. if even/odd parity
!   nlist = length of problist, hamlist, and jlist; int(Jmax-Jmin)+1
!
! OUTPUT:
!   problist(2,:) = fraction of original HF state in the norm for each J
!
!  INTERNAL:
!	TSDNmat = temporary Slater determinant Norm matrix (holds values from N_J_MK)
!	npair = number of parity states (i.e. 1 or 2) [via phf_vals module]
!
subroutine extractNorm(numOfJ,isOdd,N_J_MK,PN_J_MK,ParityTest,nlist,problist)

	use phf_vals
	use jaggedArrayType

	implicit none

!........... INPUT........................
	integer (kind=4), intent(in) :: numOfJ
	logical, intent(in) :: isOdd
	type (jaggedArray), allocatable, intent(in) :: N_J_MK(:), PN_J_MK(:)
	logical, intent(in) :: ParityTest
	integer (kind=4), intent(in) :: nlist

!............OUTPUT......................
	real (kind=8), intent(out) :: problist(2,nlist)

!............INTERNAL....................
	integer (kind=4) :: intJ, kk
	complex(kind=8), allocatable, dimension(:,:) :: TSDNmat
	integer (kind=4) :: parityOUT

	problist = 0.0d0
	do intJ = 1, numOfJ
		allocate(TSDNmat(ubound(N_J_MK(intJ)%MK,1),ubound(N_J_MK(intJ)%MK,2)))

		do parityOUT = 1, npair
			TSDNmat = cmplx(0.0d0,0.0d0)
			if (ParityTest) then
				TSDNmat = N_J_MK(intJ)%MK + (-1.0d0)**(parityOUT+1)*PN_J_MK(intJ)%MK
			else
				TSDNmat = N_J_MK(intJ)%MK
			end if

			do kk = 1, ubound(TSDNmat,1)
				problist(parityOUT,intJ) = problist(parityOUT,intJ) + dble(TSDNmat(kk,kk))
			end do
		end do

		deallocate(TSDNmat)
	end do

return; end subroutine extractNorm

!==============================================================================
!
! prints formatted problist to screen or file
!
!  INPUT:
!	numOfJ = # of beta values -> int(Jmax + 1.)
!   isOdd = logical flag for odd A nuclides
!   ParityTest = logical, .true. if even/odd parity
!   nlist = length of problist, hamlist, and jlist; int(Jmax-Jmin)+1
!   problist(2,:) = fraction of original HF state in the norm for each J
!   printTo = target of write() statements
!               6 == to screen (standard output)
!               anything >=10 is an external file  
!
!  INTERNAL:
!	jlist(:) = list of floating point J values (i.e. 0.0, 1.0, ...)
!   probsum = sum of all elements in problist(2,:) should be 1.0 (or 2.0 if ParityTest)
!
subroutine printNorm(numOfJ,isOdd,ParityTest,nlist,problist,printTo)

	implicit none

!........... INPUT........................
	integer (kind=4), intent(in) :: numOfJ
	logical, intent(in) :: isOdd
	logical, intent(in) :: ParityTest
	integer (kind=4), intent(in) :: nlist
	real (kind=8), intent(out) :: problist(2,nlist)
    integer (kind=4) :: printTo

!............INTERNAL....................
	integer (kind=4) :: intJ
    real(kind=8) :: jlist(nlist)
    real(kind=8) :: probsum

	do intJ = 1, numOfJ ! build jlist
		jlist(intJ) = dble(intJ) - 1.0d0
		if (isOdd) jlist(intJ) = jlist(intJ) + 0.5d0
	end do

    probsum = 0.d0 
    write(printTo,*)' '
    write(printTo,*)' Fraction in original HF state: Norm'
    write(printTo,*)' '
    if(ParityTest)then
	    write(printTo,*)' J    frac(+)   frac(-)'
	    write(printTo,*)'-----------------------'
	    do intJ = 1,numOfJ
		    write(printTo,2001)jlist(intJ),problist(1,intJ),problist(2,intJ)
		    probsum = probsum+problist(1,intJ)+problist(2,intJ)
	    end do
	    2001 format(f4.1,2f10.6)
    else
	    write(printTo,*)' J    frac'
	    write(printTo,*)'-----------------------'
	    do intJ = 1,numOfJ
		    write(printTo,2001)jlist(intJ),problist(1,intJ)
		    probsum = probsum+problist(1,intJ)
	    end do
    end if
    write(printTo,*)' Total of sum of trace(N) = ',probsum

return; end subroutine printNorm

