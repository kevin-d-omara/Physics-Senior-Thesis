CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  package MONTY_DENSLIB.f  
C

      subroutine makerhoij(it,np,sdf,sdi,ovlp,rhoij)
C
C  INPUT:
C    it : species = 1 proton, 2 neutron
C    np : = num of particles for this
C    sdf, sdi : final (left) and initial (right) SDs
C
C  OUTPUT:
C    ovlp: = det  sdf^+ sdi
C    rhoij : = density matrix
C           rho = sdi *( sdf^+ sdi )^-1 *sdf^+
C

      use spstate
      use system_parameters

      implicit none
      integer it
      integer np
      complex(kind = 8) :: sdf(nsps,np),sdi(nsps,np)
      complex(kind = 8) :: rhoij(nsps,nsps)
      complex(kind = 8) :: ovlp

C---------------- INTERMEDIATES ----------------

      complex(kind = 8), allocatable :: ovrtmp(:,:),X(:,:)

      complex(kind = 8) :: zzero,tmp

      integer :: i,j,k,l
C------------ FOR DECOMPOSITION

      integer,allocatable :: ipiv(:)
      integer :: info

C----------- COMPUTE  sdf^+ * sdi = ovrtmp


      zzero = (0.d0,0.d0)
      if(np == 0)then
          ovlp = (1.d0,0.d0)
		  rhoij=(0.d0,0.d0)
          return
      endif

      if(allocated(ovrtmp))deallocate(ovrtmp) !this might be inefficient; fix with pointers
      allocate(ovrtmp(np,np))

      do i = 1,np
         do j = 1,np
            tmp = zzero
            do k = 1,nsps
              tmp = tmp + dconjg(sdf(k,i))*sdi(k,j)

            enddo ! loop over k
            ovrtmp(i,j) = tmp
         enddo  ! loop over j
      enddo  ! loop over i


C------------- INVERT ovrtmp and find determinant = overlap
C              best to use LAPACK routines
c============ First LU decomposition
      if(allocated(ipiv))deallocate(ipiv)
      allocate(ipiv(np))

      call zgetrf(np,np,ovrtmp,np,ipiv,info)

      if(info /= 0)then
        print*,' problem with LU decomposition ',info
        stop
      endif

C---------- COMPUTE OVERLAP --------
      ovlp = (1.d0,0.d0)
      do i = 1,np
        ovlp = ovlp*ovrtmp(i,i)
        if(ipiv(i) /= i)ovlp = -ovlp
      enddo

C----------- INVERT. The most efficient way is to solve 
C    OVRTEMP*X = SDF^+ so that X = (OVRTMP)^-1 * SDF^+

      if(allocated(X))deallocate(x)
      allocate(X(np,nsps))
      do i = 1,np
         do j = 1,nsps
            x(i,j) = dconjg(sdf(j,i))
         enddo
      enddo
C      write(19,*)x
      call ZGETRS( 'N', np, nsps, ovrtmp, np,IPIV,X,
     &     np, INFO )
C      write(19,*)' '
C      write(19,*)x

C-------------- CONSTRUCT density matrix

      do i = 1,nsps
         do j = 1,nsps
            tmp = zzero
            do k = 1,np
                  tmp = tmp  + sdi(j,k)*X(k,i)

            enddo    ! loop over k
            rhoij(i,j) = tmp
         enddo  ! loop over j
C         write(17,101)(real(rhoij(i,j)),j=1,nsps)
101   format(6f8.4)
      enddo  ! loop over i

     

      return
      end subroutine makerhoij
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc