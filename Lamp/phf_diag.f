CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rot_and_diag
C
C  rotates and diagonalizes
C
      use system_parameters

      use spstate
      implicit none

      type wfnlist
        complex(kind = 8), allocatable :: psd(:,:),nsd(:,:)
      end type wfnlist

      type (wfnlist), allocatable :: psilist(:)

      complex(kind = 8),allocatable :: psdr(:,:),nsdr(:,:)
      complex(kind = 8),allocatable :: psdl(:,:),nsdl(:,:)
      complex(kind = 8),allocatable :: psd0(:,:),nsd0(:,:)
      complex(kind = 8),allocatable :: rhopij(:,:),rhonij(:,:)
      complex(kind = 8) :: ovlpp,ovlpn

      complex(kind = 8) :: vme

      real :: a  ! angle of rotation
      complex(kind = 8),allocatable :: drot(:,:)
      integer ia,ja

C--------- FOR TESTING ------

      integer na,amax
      integer ndim
      integer nf  ! final dimension 
      complex(kind =8),allocatable :: norm(:,:)
      complex(kind =8),allocatable :: L(:,:),Linv(:,:)
      complex(kind =8),allocatable :: hamm(:,:),h(:,:)
      integer, allocatable :: skip(:)
      real :: tol
C---------------- FOR HAMILTONIAN DIAG

      integer :: info,lwork
      real (kind = 8), allocatable :: e(:)
      real (kind = 8),allocatable :: rwork(:)
      complex(kind = 8), allocatable :: work(:)

C----------------- SET UP ARRAYS -------------------

      tol = 0.02

      print*,' '
      print*,' Enter max rotation angle in degrees '
      read*,amax
      amax = amax*3.141569/180.  ! convert to radians

      print*,' Enter # of points '
      read*,na

      ndim = na*na
 
      print*,' starting dimension= ',ndim
      allocate(psilist(ndim))
      allocate(h(ndim,ndim),norm(ndim,ndim))
      allocate(l(ndim,ndim),linv(ndim,ndim))
      allocate(skip(ndim))
      allocate(hamm(ndim,ndim))
      do ia = 1,ndim
         allocate(psilist(ia)%psd(nsps,nsps))
         allocate(psilist(ia)%nsd(nsps,nsps))
      enddo


      allocate(psdr(nsps,numprot),nsdr(nsps,numneut))
      allocate(psdl(nsps,numprot),nsdl(nsps,numneut))
      allocate(psd0(nsps,numprot),nsd0(nsps,numneut))


      call GetSherpaSD(psdr,nsdr)

      allocate(rhopij(nsps,nsps),rhonij(nsps,nsps))
      call makerhoij(1,numprot,psdr,psdr,ovlpp,rhopij)
      call makerhoij(2,numneut,nsdr,nsdr,ovlpn,rhonij)

C      print*,' Overlaps ',ovlpp,ovlpn

      call TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
      print*,' starting  < H > = ',vme

C----------- ROTATE AND STORE

      allocate(drot(nsps,nsps))

      do ia = 1,na
      a = (ia-1)*amax/(na-1)
      call makeZrot(nsps,a,drot)
      call rotate(nsps,numprot,drot,psdr,psd0)
      call rotate(nsps,numneut,drot,nsdr,nsd0)
        do ja = 1,na
           a = (ja-1)*amax/(na -1)
          call makeYrot(nsps,a,drot)
          call rotate(nsps,numprot,drot,psd0,psdl)
          call rotate(nsps,numneut,drot,nsd0,nsdl)
          psilist(ja+(ia-1)*na)%psd = psdl
          psilist(ja+(ia-1)*na)%nsd = nsdl
        enddo
      enddo
C-------------- NOW MAKE MATRICES 

      do ia = 1,ndim
          psdr = psilist(ia)%psd
          nsdr = psilist(ia)%nsd
          do ja = ia,ndim
             psdl = psilist(ja)%psd
             nsdl = psilist(ja)%nsd

             call makerhoij(1,numprot,psdl,psdr,ovlpp,rhopij)
             call makerhoij(2,numneut,nsdl,nsdr,ovlpn,rhonij)
             norm(ia,ja) = ovlpp*ovlpn
             norm(ja,ia) = dconjg(norm(ia,ja))
             call TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
             hamm(ia,ja) = vme
             hamm(ja,ia) = dconjg(vme)

          enddo
      enddo

      call choleskymaster(ndim,ndim,norm,tol,.true.,L,Linv,skip)

      call choleskyinvert(ndim,ndim,nf,hamm,h,Linv,skip)

C--------  SET UP FOR DIAGONALIZATION 

      allocate(e(nf))
      lwork = 2*nf - 1
      allocate(work(lwork))
      allocate(rwork(3*nf -2  ) )
      call zheev('n','u',nf,h,ndim,e,work,lwork,rwork,info)

      print*,e(1)
      return

      end subroutine rot_and_diag


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC