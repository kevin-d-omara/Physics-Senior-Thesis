	module phf_vals
C
C  projected hartree-fock norm and operator matrix values
C
	implicit none

C------------ CREATE A DEFINED TYPE----------------

	type phf_type
		real (kind = 8) :: j
		integer :: tdim
		complex (kind = 8), allocatable, dimension(:,:) :: h
          complex (kind = 8), allocatable, dimension(:,:) :: n
	end type phf_type
	
C
C    PHF - no parity, PHFP - raw parity, PHFM - Had an idea for this...now its gone
C
	  type (phf_type), allocatable :: phf(:), phfp(:), phfm(:)
C
C    JLOC - ang. mom location index (e.g. J = 0.5, jloc = 0, J=1.5, jloc = 1)
C    used for the allocated array of phf_type index
C    NPAIR - used to allocate parity/no-parity e-value/e-vec arrays
C	         NPAIR = 1 (no parity), NPAIR = 2 (parity)
C
       integer :: jloc,npair
	
	end module phf_vals
	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	module phf_allocation
	
	contains
C
C    MEMORY ALLOCATION FOR PHF STRUCTURES
C	
	subroutine phf_alloc(jmin,jmax,pairtest)
		
	     use phf_vals
		
		implicit none
		
		real(kind=8),intent(in) :: jmin,jmax
		logical, intent(in) :: pairtest
		integer :: numj,i,intj
		real :: tempj
		
		numj = int(jmax - jmin)
				
		allocate(phf(0:numj))
		if (pairtest) allocate(phfp(0:numj))
		if (pairtest) allocate(phfm(0:numj))
		!	Flag for systematic deallocation

		phf(:)%tdim = numj
		if (pairtest) phfp(:)%tdim = numj
		if (pairtest) phfm(:)%tdim = numj
		
		tempj = jmin
		do i = 0, numj
			intj = int(2.0d0*tempj) + 1
			allocate(phf(i)%h(intj,intj),phf(i)%n(intj,intj))

               if (pairtest) allocate(phfp(i)%h(intj,intj))
               if (pairtest) allocate(phfp(i)%n(intj,intj))
               
               if (pairtest) allocate(phfm(i)%h(intj,intj))
               if (pairtest) allocate(phfm(i)%n(intj,intj))
			
               phf(i)%j = tempj
			if (pairtest) phfp(i)%j = tempj
			if (pairtest) phfm(i)%j = tempj
			tempj = tempj + 1.0d0
		end do		
	
	end subroutine
C
C    MEMORY DEALLOCATION FOR PHF STRUCTURES
C
	
	subroutine phf_dealloc
	
		use phf_vals
	
		implicit none

		integer :: numj,i
	
		if (allocated(phf)) then
			numj = phf(0)%tdim
			do i = 0, numj
               			deallocate(phf(i)%h,phf(i)%n)
			end do
               deallocate(phf)
		end if

          if (allocated(phfp)) then
               numj = phfp(0)%tdim
               do i = 1, numj
                    deallocate(phfp(i)%h,phf(i)%n)
               end do
               deallocate(phfp)
          end if

          if (allocated(phfm)) then
               numj = phfm(0)%tdim
               do i = 1, numj
                    deallocate(phfm(i)%h,phfm(i)%n)
               end do
               deallocate(phfm)
          end if

	end subroutine phf_dealloc
	
	end module phf_allocation

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module sporbit
C
C  single-particle ORBIT information
C
      implicit none

      integer numorb          ! # of s.p. orbits
      integer numruns

C------------ CREATE A DEFINED TYPE----------------
      type orb
        integer :: nr       ! radial quantum number
        integer :: j        ! 2 x j
        integer :: l        ! L
        integer :: w        ! excitation 
        integer :: spstart  ! where these orbits correspond to 
                            ! start in spqn

      end type orb

      type (orb),allocatable :: orbqn(:)

      end module sporbit

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module spstate
C
C  single-particle STATE information
C
      implicit none

      integer nsps          ! # of s.p. states

C------------ CREATE A DEFINED TYPE----------------
      type spst
        integer :: nr       ! radial quantum number
        integer :: j        ! 2 x j
        integer :: m        ! 2 x jz
        integer :: l        ! L
        integer :: w        ! excitation 
        integer :: par      ! parity
        integer :: orb      ! orbit label
      end type spst

      type (spst),allocatable :: spsqn(:)

      end module spstate

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module pairdef

      type pair_qn
        integer :: M
        integer :: par
        integer :: W
        integer :: indx
        integer :: ia,ib
      end type pair_qn

      end module pairdef
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module interaction

      use pairdef
      real,allocatable :: spe(:)


      integer :: nvdim  ! dimension

      integer :: ntbme     ! # of two-body matrix elements

      integer :: norb_allow  ! # of allowed orbits
      integer,allocatable :: orblist(:)  

C----------- coupled  TBMES -----------------
      type vjts
        integer  :: Jmin,Jmax
        real,allocatable :: v(:,:)
        integer  :: orb(4)

      end type vjts

      type (vjts),allocatable :: vtbme(:)  

C---------- uncoupled TBMEs ------------------


      integer :: nmatpp,nmatnn
      integer :: nmatpn
      integer :: nmatXX

      real,allocatable  :: hmatxx(:,:)
      integer, allocatable :: hmatorbxx(:,:,:)


      real, allocatable :: hmatpn(:)
      integer, allocatable :: hmatorbPN(:,:)

C      real,allocatable  :: hdiag(:)  ! diagonal matrix elements

C------------ coding of the pairs -------------
C  take states i,j,k,l  from list in hspsqn
C
C  We encode a_i a_j, or a^+_i a^+_j, as follows:
C
C  assume i >= j
C
C  then  (i,j) =>  i(i-1)/2 +j 
C

C-----------------------------

C------ uncoupled SPEs

      real, allocatable :: speunX(:,:,:),speunXshift(:,:,:)

C------------ SPE shifts

      real, allocatable :: speshift(:)  !ULTIMATELY ALLOW FOR different for p,n

      end module interaction
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module system_parameters

      implicit none

      integer numprot,numneut,numsd
      integer jz
      integer iparity
      integer             :: Nrank  ! rank of interaction, e.g. 2-particle, 3-particle

      end module system_parameters

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module decomposition
      implicit none

C------------max K for rho(i,j)      
      integer :: kmax

C----------- list of max_k for various rho(i,j)

      integer,allocatable :: maxk(:,:), mink(:,:)

C----------- list of rho(i,j) as a function of K
      type rhox
         integer :: nrhos  ! # of density operators
         integer, allocatable :: i(:),j(:)
      end type rhox      
      type (rhox), allocatable :: rholist(:)  ! function of K

C---------- basis for PQ operators

      type protoPQ

         complex (kind =8), allocatable :: P(:,:)  ! function of i,j
         complex (kind =8), allocatable :: Q(:,:)  ! function of i,j

      end type protoPQ

C---------- pandya-transformed hamiltonian 
      type pandyamat
         real, allocatable :: h(:,:)
         real, allocatable :: lam(:)
         real, allocatable :: vec(:,:)
         complex,allocatable :: s(:)  ! sign function, equals 1 or i
         integer :: nops
         type (protoPQ), allocatable :: op(:,:)  ! function of nops, M >=0
      end type pandyamat
      type (pandyamat), allocatable :: E_k(:)  !function of K



      end module decomposition

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module auxfields

      implicit none

      type fields
         complex(kind=8),allocatable :: sig(:,:)  ! of alpha,m
         complex(kind=8),allocatable :: tau(:,:)  ! of alpha,m

      end type fields

      type(fields),allocatable :: aux_k(:)  ! of K

      type(fields),allocatable :: shift_k(:)  ! of K! 

C........... TIME-DEPENDENT FIELDS .........................     

      type stfields
          type(fields),allocatable :: fk(:)  ! of K

      end type stfields

      type(stfields), pointer :: aux_t(:)    ! of it = time

      end module auxfields


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module psis
C
C  time-dependent wavefunctions
C
      implicit none


      type sdmaster
          complex(kind=8),pointer :: psd(:,:),nsd(:,:)
      end type sdmaster

      type(sdmaster), pointer :: psir(:),psil(:)  ! of it = time


      end module psis

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module timestep

      implicit none
      real :: beta  ! total beta

      integer :: nt ! # time steps
      real    :: dbeta

      integer npade   ! order for pade approximant

      end module timestep

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module errortests

      implicit none
      
      logical :: diagp
      character(LEN = 100) :: path

      end module errortests

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module jaggedArrayType

      type jaggedArray
          complex(kind=8), allocatable, dimension(:,:) :: MK
      end type jaggedArray

      end module jaggedArrayType
