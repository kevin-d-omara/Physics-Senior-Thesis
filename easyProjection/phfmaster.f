CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      program montymaster
      
      use errortests 

      implicit none
      
      path = "undefined"
      path = TRIM(path)

      call J_Main

      end

C
C  master calling program for Auxiliary Field Monte Carlo
C  CWJ / SDSU / September 2007
C 
      subroutine testset

      use system_parameters

      call get_sp_info

      call setup4tbmes
      call readvtbme
      call uncoupleXXmaster
      call undoSPE
      call uncouplePNmaster
      print*,' Enter Z, N '
      read*,numprot,numneut

      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine setupwfns
C
C  subroutine to test reading in (a) slater determinant(s)
C  and computing the matrix element
C

      use system_parameters
      use spstate
      use psis
      implicit none

      integer it

      complex(kind = 8),allocatable :: psd(:,:),nsd(:,:)

      complex(kind = 8),allocatable :: rhopij(:,:),rhonij(:,:)
      complex(kind = 8) :: ovlpp,ovlpn

      complex(kind = 8) :: vme

      allocate(psd(nsps,numprot),nsd(nsps,numneut))
      call GetSherpaSD(psd,nsd)

      allocate(rhopij(nsps,nsps),rhonij(nsps,nsps))
      call makerhoij(1,numprot,psd,psd,ovlpp,rhopij)
      call makerhoij(2,numneut,nsd,nsd,ovlpn,rhonij)
      print*,' Overlaps ',ovlpp,ovlpn
C      write(17,*)rhopij

      call TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
      print*,' < H > = ',vme

C-------------- NOW STORE

      !psir(0)%psd = psd
      !psil(0)%psd = psd
      !psir(0)%nsd = nsd
      !psil(0)%nsd = nsd

      return
      end subroutine setupwfns

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine TestReadAndOverlap
C
C  subroutine to test reading in (a) slater determinant(s)
C  and computing overlap(s)
C

      use system_parameters
      use spstate
      implicit none

      complex(kind = 8),allocatable :: psdi(:,:),psdf(:,:)
      complex(kind = 8),allocatable :: nsdi(:,:),nsdf(:,:)

      complex(kind = 8),allocatable :: rhopij(:,:),rhonij(:,:)
      complex(kind = 8) :: ovlpp,ovlpn

      allocate(psdi(nsps,numprot),nsdi(nsps,numneut))
      call GetSherpaSD(psdi,nsdi)

      allocate(rhopij(nsps,nsps),rhonij(nsps,nsps))
      call makerhoij(1,numprot,psdi,psdi,ovlpp,rhopij)
      call makerhoij(2,numneut,nsdi,nsdi,ovlpn,rhonij)

      print*,' Overlaps ',ovlpp,ovlpn

      allocate(psdf(nsps,numprot),nsdf(nsps,numneut))
      call GetSherpaSD(psdf,nsdf)

      call makerhoij(1,numprot,psdf,psdi,ovlpp,rhopij)
      call makerhoij(2,numneut,nsdf,nsdi,ovlpn,rhonij)

      print*,' Overlaps ',ovlpp,ovlpn


      return
      end subroutine TestReadAndOverlap

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine TestReadAndHam
C
C  subroutine to test reading in (a) slater determinant(s)
C  and computing the matrix element
C

      use system_parameters
      use spstate
      implicit none

      complex(kind = 8),allocatable :: psd(:,:),nsd(:,:)

      complex(kind = 8),allocatable :: rhopij(:,:),rhonij(:,:)
      complex(kind = 8) :: ovlpp,ovlpn

      complex(kind = 8) :: vme

      allocate(psd(nsps,numprot),nsd(nsps,numneut))
      call GetSherpaSD(psd,nsd)

      allocate(rhopij(nsps,nsps),rhonij(nsps,nsps))
      call makerhoij(1,numprot,psd,psd,ovlpp,rhopij)
      call makerhoij(2,numneut,nsd,nsd,ovlpn,rhonij)

      print*,' Overlaps ',ovlpp,ovlpn
C      write(17,*)rhopij

      call TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
      print*,' < H > = ',vme

      return
      end subroutine TestReadAndHam

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 