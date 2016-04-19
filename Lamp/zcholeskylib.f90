!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  zcholeskylib
!
!  routines for double precision complex cholesky decomposition
!
!  including deletion of the basis
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine testcholesky
!
!  master subroutine to test rotations and cholesky
!
      use system_parameters

      use spstate
      implicit none

      type wfnlist
        complex(kind = 8), allocatable :: psd(:,:),nsd(:,:)
      end type wfnlist

      type (wfnlist), allocatable :: psilist(:)

      complex(kind = 8),allocatable :: psdr(:,:),nsdr(:,:)
      complex(kind = 8),allocatable :: psdl(:,:),nsdl(:,:)

      complex(kind = 8),allocatable :: rhopij(:,:),rhonij(:,:)
      complex(kind = 8) :: ovlpp,ovlpn

      complex(kind = 8) :: vme

      real :: a  ! angle of rotation
      complex(kind = 8),allocatable :: drotZ(:,:)
      integer ia,ja

!--------- FOR TESTING ------

      integer na
      complex(kind =8),allocatable :: h(:,:),norm(:,:)
      complex(kind =8),allocatable :: L(:,:),Linv(:,:)
      complex(kind =8),allocatable :: out(:,:)
      integer, allocatable :: skip(:)
      real :: tol
	  logical :: posdef

!----------------- SET UP ARRAYS -------------------

      tol = 0.1

      print*,' Enter # of points '
      read*,na
      allocate(psilist(na))
      allocate(h(na,na),norm(na,na),l(na,na),linv(na,na))
      allocate(skip(na))
      allocate(out(na,na))
      do ia = 1,na
         allocate(psilist(ia)%psd(nsps,nsps))
         allocate(psilist(ia)%nsd(nsps,nsps))
      enddo


      allocate(psdr(nsps,numprot),nsdr(nsps,numneut))
      allocate(psdl(nsps,numprot),nsdl(nsps,numneut))

      call GetSherpaSD(psdr,nsdr)

      allocate(rhopij(nsps,nsps),rhonij(nsps,nsps))
      call makerhoij(1,numprot,psdr,psdr,ovlpp,rhopij)
      call makerhoij(2,numneut,nsdr,nsdr,ovlpn,rhonij)

      print*,' Overlaps ',ovlpp,ovlpn

      call TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
      print*,' < H > = ',vme

!----------- ROTATE AND STORE

      allocate(drotZ(nsps,nsps))

      psilist(1)%psd = psdr
      psilist(1)%nsd = nsdr
      do ia = 2,na
      a = (ia-1)*90./na
      call makeZrot(nsps,a,drotZ)
      call rotate(nsps,numprot,drotZ,psdr,psdl)
      call rotate(nsps,numneut,drotZ,nsdr,nsdl)
      psilist(ia)%psd = psdl
      psilist(ia)%nsd = nsdl

      enddo
!-------------- NOW MAKE MATRICES 

      do ia = 1,na
          psdr = psilist(ia)%psd
          nsdr = psilist(ia)%nsd
          do ja = ia,na
             psdl = psilist(ja)%psd
             nsdl = psilist(ja)%nsd

             call makerhoij(1,numprot,psdl,psdr,ovlpp,rhopij)
             call makerhoij(2,numneut,nsdl,nsdr,ovlpn,rhonij)
             norm(ia,ja) = ovlpp*ovlpn
             norm(ja,ia) = dconjg(norm(ia,ja))
          enddo
      enddo

      call choleskymaster(na,na,norm,tol,.true.,L,Linv,skip,posdef)
!      call choleskyinvert(na,na,norm,out,Linv,skip)
      do ia = 1,na
         if(skip(ia) ==1)print*,ia,out(ia,ia),out(1,ia)
      enddo

      return

      end subroutine testcholesky


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine choleskyinvert(np,n,nf,A,B,Linv,skip)
!
!  applies Linv to A to get B
!  B = Linv A Linv^dagger
!
!  also eliminates "skipped" rows
!
      implicit none
      integer NP,N 	! dimension of all array
      integer nf        ! final dimension after skipping
      complex(kind = 8) :: A(NP,NP),B(Np,NP),Linv(NP,np)
      integer ::           skip(NP)
      real    ::           tol
      logical ::           invflag

      

      integer i,j,k,m
      integer icount,jcount

      complex(kind =8) :: zsum,ztmp

      nf = 0
      do i = 1,n
        nf = nf+skip(i)
      enddo

      b = (0.d0,0.d0)
      jcount = 0
      do j = 1,N
        if(skip(j) == 0)cycle
        jcount = jcount + 1
        icount = 0
        do i = 1,j
        if(skip(i)==0)cycle
        icount = icount + 1
        zsum = (0.d0,0.d0)
        do m = 1,j
           ztmp = dconjg(linv(j,m))*skip(m)
           do k = 1,i
             zsum = zsum + linv(i,k)*a(k,m)*ztmp*skip(k)
           enddo
        enddo
        b(icount,jcount) = zsum
        b(jcount,icount) = dconjg(zsum)
        enddo
      enddo

      return
      end subroutine choleskyinvert
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine choleskymaster(np,n,A,tol,invflag,L, Linv, skip,posdef)
!
!  master subroutine for Cholesky
!
! INPUT
!    NP = declared dimension of all arrays
!    N =  actual dimension ( usually = NP)
!    A = original (complex*8) matrix to be decomposed
!    tol = tolerance for diagonal (real*4)
!    invflag = .true. to compute inverse
! OUTPUT
!    L    = lower triangular cholesky (complex*8)
!    Linv = lower triangular of inverse
!    skip = integer array keeps track of what to keep
!
!  SUBROUTINES CALLED:
!
!  zcholeskystep
!
      implicit none
	  interface
	      subroutine zcholeskystep(n,i,A,tol,invflag,L,Linv,skip,posdef)  ! INTERFACE
	      use errortests
	      implicit none
	      integer N 	! dimension of all array
	      complex(kind = 8) :: A(N,N), L(N,N),Linv(N,n)
	      integer i         ! which step we're at
	      integer ::           skip(N)
	      real    ::           tol
	      logical ::           invflag,posdef
	      end subroutine zcholeskystep ! INTERFACE
	  end interface
	  
      integer NP,N 	! dimension of all array
      complex(kind = 8) :: A(NP,NP), L(NP,NP),Linv(NP,np)
      integer ::           skip(NP)
      real    ::           tol
      logical ::           invflag, posdef

      integer i
      integer nsave

      complex(kind = 8),allocatable :: b(:,:)
      integer j,k

      posdef = .true.

      do i = 1,n
        call zcholeskystep(np,i,A,tol,invflag,L,Linv,skip,posdef)
	if (posdef .EQV. .false.) return
      enddo

      nsave = 0
      do i = 1,n
        nsave = nsave+skip(i)
      enddo
      return

      end subroutine choleskymaster

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine zcholeskystep(n,i,A,tol,invflag,L,Linv,skip,posdef)

      use errortests

!
!    COMPUTES CHOLESKY DECOMPOSITION FOR DOUBLE COMPLEX
!    ALSO CAN COMPUTE INVERSE SIMULTANEOUSLY
!    ALLOWS FOR (NEARLY) LINEARLY DEPENDENT BASIS BY 
!    "THROWING OUT" ROWS/COLUMNS VIA skip ARRAY
!
!    WRITTEN SO THAT THIS CAN BE DONE ONE ROW/COLUMN AT A TIME
!    ALLOWS TO FOLLOW CONVERGENCE, OR, ALTERNATELY,
!    ENLARGE SIZE OF YOUR SPACE WITHOUT REDOING WORK
!
!  CWJ /SDSU / Feb 2008
!
! INPUT
!    N = declared dimension of all arrays
!    i = current step  <= N
!    A = original (complex*8) matrix to be decomposed
!    tol = tolerance for diagonal (real*4)
!    invflag = .true. to compute inverse
! OUTPUT
!    L    = lower triangular cholesky (complex*8)
!    Linv = lower triangular of inverse
!    skip = integer array keeps track of what to keep
!
!   SUBROUTINES CALLED:
!

      implicit none
      integer N 	! dimension of all array
      complex(kind = 8) :: A(N,N), L(N,N),Linv(N,n)
      integer i         ! which step we're at
      integer ::           skip(N)
      real    ::           tol
      logical ::           invflag,posdef

!------------------ USED INTERNAL

      integer j,k
      complex(kind =8) ::  sum
      real   (kind =8) ::  rsum

!-------------- ERROR TRAP FOR DIAGONAL ------------
      diagp = .true.

      if(dreal(A(i,i)) < -tol )then
         !print*,i,' real part of diagonal < 0 ',A(i,i)
         diagp = .false.
	 return
      endif
      if(abs(dimag(A(i,i))) > tol )then
         !print*,i,' imag part of diagonal not 0 ',A(i,i)
         diagp = .false.
         return
      endif

!--------------END ERROR TRAP ----------------

!-------------- COMPUTE OFF-DIAGONALS L(i,j), j < i

      if(i > 1)then
        do j = 1,i-1
          if(skip(j) == 0)cycle
          sum = A(i,j)
          do k = 1,j-1
             sum = sum - L(i,k)*dconjg(L(j,k))*skip(k)
          enddo
          L(i,j) = sum/L(j,j)
        enddo  
      endif

!------------- CHECK TO SKIP -----------------

      rsum = dreal(A(i,i))
      if(i > 1)then
        do j = 1,i-1
           rsum = rsum - abs(L(i,j))*abs(L(i,j))*skip(j)
        enddo
      endif
      if(rsum < -tol)then
           !print*,' Matrix not pos def '
           !print*, i,rsum
	   posdef = .false.
           return
      endif

      if(rsum < tol)then
        skip(i) = 0
        return
      else
        skip(i) = 1
      endif

      L(i,i) = cmplx(dsqrt(rsum),0.d0)

      if(.not.invflag)return
!----------------- COMPUTE INVERSE -------------

      Linv(i,i) = 1.d0 / L(i,i)
      if(i == 1) return

      do j = i-1,1,-1
        if(skip(j) == 0)cycle
        sum = 0.d0
        do k = j+1,i
           sum = sum + Linv(i,k)*L(k,j)*skip(k)
        enddo
        Linv(i,j) = -sum/L(j,j)
      enddo

      return
      end subroutine zcholeskystep