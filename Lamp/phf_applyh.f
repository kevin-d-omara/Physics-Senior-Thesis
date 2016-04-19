!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  package monty_applyh
!C
!C  routines to compute matrix elements
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)

      use system_parameters
      implicit none
      integer nsps

      complex(kind = 8) :: rhopij(nsps,nsps),rhonij(nsps,nsps)
      complex(kind = 8) :: ovlpp,ovlpn

      complex(kind = 8) :: vme, vtbmePP,vtbmeNN,vtbmePN

      
      vtbmePP = (0.d0,0.d0)
      vtbmeNN = (0.d0,0.d0)
      vtbmePN = (0.d0,0.d0)



      if(numprot > 1)then
         call TBMExx(1,nsps,rhopij,vtbmePP)
      endif
      if(numneut > 1 ) then
         call TBMExx(1,nsps,rhonij,vtbmeNN)

      endif

      if(numprot > 0)then
         call applySPE(1,nsps,rhopij,vtbmePP)
      endif
      if(numneut > 0)then
         call applySPE(2,nsps,rhonij,vtbmeNN)

      endif

      if(numprot > 0 .and. numneut > 0)then
           call TBMEpn(nsps,rhopij,rhonij,vtbmePN)

      endif

      vme = vtbmePP+vtbmeNN + vtbmePN

      vme = vme*ovlpp*ovlpn

      return
      end subroutine TBMEmaster

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine TBMExx(it,nsps,rhoij,vtbmeXX)
!C
!C  compute PP/NN matrix elements
!C
      use interaction
      implicit none

      integer it    ! species 1 = protons 2 = neutrons
      integer nsps
      complex(kind = 8):: rhoij(nsps,nsps)
      complex(kind = 8) :: vtbmeXX
      complex(kind = 8) :: zsum
      real(kind = 4) :: vtmp
      integer itbme
      integer a,b,c,d

!C--------------- loop over matrix elements
!$OMP PARALLEL shared(hmatorbXX,hmatXX), private(a,b,c,d,vtmp,zsum)
!$OMP  do schedule(static) reduction(+:vtbmeXX)
      do itbme = 1,nmatXX
          a = hmatorbXX(it,itbme,1)  
          b = hmatorbXX(it,itbme,2)
          c = hmatorbXX(it,itbme,3)
          d = hmatorbXX(it,itbme,4)

          vtmp = hmatXX(it,itbme)
!C  ------------- FIND ALL PERMUTATIONS 

          zsum = rhoij(a,c)*rhoij(b,d) - rhoij(a,d)*rhoij(b,c)
     &          +rhoij(c,a)*rhoij(d,b) - rhoij(d,a)*rhoij(c,b)
          vtbmeXX = vtbmeXX + vtmp*zsum
      enddo  ! loop over itbme
!$OMP END  DO 
!$OMP END PARALLEL

      return
      end subroutine TBMExx

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine TBMEpn(nsps,rhopij,rhonij,vtbmePN)
!C
!C  compute PP/NN matrix elements
!C

      use interaction
      implicit none

      integer nsps
      complex(kind = 8):: rhopij(nsps,nsps),rhonij(nsps,nsps)
      complex(kind = 8) :: vtbmePN
      complex(kind = 8) :: zsum
      real(kind = 4) :: vtmp
      integer itbme
      integer a,b,c,d

!C--------------- loop over matrix elements
!$OMP PARALLEL shared(hmatorbPN,hmatPN), private(a,b,c,d,vtmp,zsum)
!$OMP  do schedule(static) reduction(+:vtbmePN)
      do itbme = 1,nmatPN
          a = hmatorbPN(1,itbme)  
          b = hmatorbPN(2,itbme)
          c = hmatorbPN(3,itbme)
          d = hmatorbPN(4,itbme)

          vtmp = hmatPN(itbme)
C  ------------- FIND ALL PERMUTATIONS 

          zsum = rhopij(a,c)*rhonij(b,d) 
          vtbmePN = vtbmePN + vtmp*zsum
      enddo  ! loop over itbme
!$OMP END  DO 
!$OMP END PARALLEL	  
      return
      end subroutine TBMEpn
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine applySPE(it,nsps,rhoij,vmeX)
!
!  adds single-particle energies
!
      use interaction
      implicit none

      integer it    ! species 1 = protons 2 = neutrons
      integer nsps
      complex(kind = 8):: rhoij(nsps,nsps)
      complex(kind = 8) :: vmeX
      complex(kind = 8) :: zsum
      real(kind = 4) :: vtmp
      integer itbme
      integer a,b

      do a = 1,nsps
         do b = 1,nsps
           vmeX = vmeX + speunX(b,a,it)*rhoij(b,a)  ! ASSUME SPE symmetric
         enddo  ! loop over b
      enddo  ! loop over a

      return
      end subroutine applySPE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC