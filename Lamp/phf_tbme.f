CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C library Monty_TBME.f
C
C  routines to uncouple tbme in compact form
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine uncoupleXXmaster

      implicit none

C-------------- INTERFACES ---------------
      interface

      subroutine countcreateuncoupledpairs(flag,m,n,a2,a1)
        implicit none
        logical flag
        integer m,n
        integer, pointer :: a2(:,:),a1(:)
      end subroutine countcreateuncoupledpairs

      subroutine countuncoupledtbmeXX(n,a)
        implicit none
        integer n
        integer, pointer :: a(:)
      end subroutine countuncoupledtbmeXX

      subroutine untbmeXX(m,n,a2,a1)
        implicit none
        integer m,n
        integer, pointer :: a2(:,:),a1(:)
      end subroutine untbmeXX

      end interface
C-------------- END INTERFACES --------------
      integer mmax
      integer nuncpairs
      integer, pointer :: sppair(:,:),mvpair(:)
      integer nutbmeXX
      print*,' about to decouple '
      call countcreateuncoupledpairs(.true.,mmax,nuncpairs,
     & sppair,mvpair)

      call countcreateuncoupledpairs(.false.,mmax,nuncpairs,
     & sppair,mvpair)

      call countuncoupledtbmeXX(nuncpairs,mvpair)
      call untbmeXX(mmax,nuncpairs,sppair,mvpair)
      deallocate(sppair,mvpair)
      print*,' decoupled '
      return
      end subroutine uncoupleXXmaster

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine countcreateuncoupledpairs(countflag,jmax,nuncpairs,
     & sppair,mvpair)

      use spstate
      use sporbit
      use interaction

      implicit none

      integer jmax
      integer j
      integer itbme
      integer asp,bsp
      integer ma,mb
      integer nuncpairs
      real tol
      logical countflag

      integer, pointer :: sppair(:,:),mvpair(:)

      tol = 1.0e-5

C-------------------- find maxJ = max M

      if(countflag)then
      jmax = -99
      do itbme = 1,nvdim   ! loop over coupled TBMEs
        do J = vtbme(itbme)%jmin,vtbme(itbme)%jmax
           if( abs(vtbme(itbme)%v(j,1)) > tol)jmax = max(j,jmax)
        enddo  ! loop over J
      enddo  ! loop over itbme

C------------------- 
C      if(jmax == -99)then
C         stop  ! need to figure this out
C      endif
      endif
C........... COUNT UP # OF UNCOUPLED PAIRS 

      nuncpairs = 0
      do asp = 2,nsps 
        ma = spsqn(asp)%m
        do bsp = 1,asp-1
            mb = spsqn(bsp)%m
            if(abs(ma+mb) <= 2*jmax)then
               nuncpairs = nuncpairs + 1
               if(.not.countflag)then
C                  print*,nuncpairs,asp,bsp
                  sppair(nuncpairs,1) = asp
                  sppair(nuncpairs,2) = bsp
                  mvpair(nuncpairs) = (ma+mb)/2
               endif
            endif
        enddo  ! loop over bsp

      enddo  ! loop over asp
      print*,' there are ',nuncpairs,' uncoupled pairs '
      if(countflag)then
            allocate(sppair(nuncpairs,2),mvpair(nuncpairs))
      endif
      return
      end subroutine countcreateuncoupledpairs

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine countuncoupledtbmeXX(nuncpairs,mvpair)
      
      use spstate
      use interaction

      implicit none

      integer nuncpairs

      integer m
      integer i,j
      integer,pointer :: mvpair(:)

      nmatXX = 0

      do i = 1,nuncpairs
        m = mvpair(i)
        do j = i,nuncpairs
          if(m == mvpair(j))nmatXX = nmatXX +1
        enddo
      enddo  

      print*,' There are ',nmatxx,' PP/NN uncoupled TBMEs '
      if (allocated(hmatxx).or.allocated(hmatorbxx)) then
            deallocate(hmatxx,hmatorbxx)
      end if
      allocate(hmatxx(1,nmatxx),hmatorbxx(1,nmatxx,4))
      return
      end subroutine countuncoupledtbmeXX

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine untbmeXX(jmax,nuncpairs,sppair,mvpair)

      use interaction
      use spstate

      implicit none
      integer jmax
      integer nuncpairs
      integer,pointer :: sppair(:,:),mvpair(:)
      
      integer a,b,c,d
      integer asp,bsp,csp,dsp
      integer ma,mb,mc,md
      integer ja,jb,jc,jd
      integer m
      
      integer imatXX
      integer i,j
      integer ipair,jpair
      integer indx
      integer JJ
      real vtmp
      real zeta,cleb

      imatXX = 0
      do i = 1,nuncpairs
         m = mvpair(i)
         asp = sppair(i,1) 
         bsp = sppair(i,2)

         a = spsqn(asp)%orb 
         b = spsqn(bsp)%orb

         if(b > a)then
           print*,' wrong order ',a,b
           stop
         endif

         ja = spsqn(asp)%j
         jb = spsqn(bsp)%j
         ma = spsqn(asp)%m
         mb = spsqn(bsp)%m
         ipair = a*(a-1)/2 + b
C         write(18,*)i,asp,bsp,a,b,spsqn(a)%m,spsqn(b)%m
         do j = 1,i
           if(m == mvpair(j))then
              imatXX = imatXX + 1

C----------------- EXTRACT THE CORRESPONDING STATES

              csp = sppair(j,1)  
              dsp = sppair(j,2)

              c = spsqn(csp)%orb
              d = spsqn(dsp)%orb
              if(d > c)then
                print*,' wrong order c d ',c,d
                stop
              endif

              jc = spsqn(csp)%j
              jd = spsqn(dsp)%j
              mc = spsqn(csp)%m
              md = spsqn(dsp)%m
              jpair = c*(c-1)/2+d
C              write(17,*)i,j,':',a,b,c,d

              if(jpair <= ipair)then 
                 indx = ipair*(ipair-1)/2+jpair
              else
                indx = jpair*(jpair-1)/2+ipair
              endif

              vtmp = 0.0
              do JJ = vtbme(indx)%jmin,min(jmax,vtbme(indx)%jmax)
                 vtmp = vtmp+ vtbme(indx)%v(jj,1) 
     &           *zeta(a,b)*zeta(c,d)
     &           *cleb(ja,ma,jb,mb,2*jj,2*m)*cleb(jc,mc,jd,md,2*jj,2*m)
              enddo  ! loop over j

C---------------- KEEP FROM DOUBLE COUNTING FOR DIAGONAL 

              if(asp == csp .and. bsp == dsp)vtmp = vtmp*.5

              hmatXX(1,imatXX) = vtmp
              hmatorbXX(1,imatxx,1) = asp
              hmatorbXX(1,imatxx,2) = bsp
              hmatorbXX(1,imatxx,3) = csp
              hmatorbXX(1,imatxx,4) = dsp

           endif

         enddo ! loop over j
      enddo    ! loop over i
      return
      end subroutine untbmeXX

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       subroutine uncouplePNmaster

       use spstate
       implicit none


       call countcreatetbmePN(.true.)

       call countcreatetbmePN(.false.)

       return
      
       end subroutine uncouplePNmaster

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

       subroutine countcreatetbmePN(countflag)

       use spstate
       use interaction

       implicit none

       logical :: countflag
C       integer,pointer ::  mXops(:,:),parXops(:,:)

       integer nct

       integer pa,pc,nb,nd
       integer a,b,c,d
       integer jpa,jpc,mpa,mpc
       integer jnb,jnd,mnb,mnd

       integer Mp, Mn
       integer parP,parN
       integer Jmin,Jmax,J,m

       logical phaseab,phasecd
       integer fact0,fact1

       integer pair1,pair2,indx

       real vtmp

C------------ FUNCTIONS CALL
       real cleb,zeta

       nct = 0
C------------- LOOP OVER PROTON OPS

       do pa = 1,nsps
          a = spsqn(pa)%orb
          jpa = spsqn(pa)%j
          mpa = spsqn(pa)%m

          do pc = 1,nsps
             c = spsqn(pc)%orb
             jpc = spsqn(pc)%j
             mpc = spsqn(pc)%m                        
             Mp = mpa - mpc
             parP = spsqn(pa)%par*spsqn(pc)%par

C------------- loop over NEUTRON OPS
             do nb = 1,nsps
                  b = spsqn(nb)%orb
                  jnb = spsqn(nb)%j
                  mnb = spsqn(nb)%m
                  do nd = 1,nsps
                    d = spsqn(nd)%orb
                    jnd =spsqn(nd)%j
                    mnd = spsqn(nd)%m
                    Mn = mnb - mnd
                    if( Mp + Mn /= 0)cycle

                    parN = spsqn(nb)%par*spsqn(nd)%par
                    if( parN /= parP)cycle

C----------------------- FIND INDEX OF MATRIX ELEMENT --------
C                        plus phases
C               standard is a >= b, c >= d
                    if(a >= b)then
                       pair1 = a*(a-1)/2 + b
                       phaseab = .false.
                    else
                       pair1 = b*(b-1)/2 + a			
                       phaseab = .true.
                    endif

                    if(c >= d)then
                       pair2 = c*(c-1)/2 + d
                       phasecd = .false.
                    else
                       pair2 = d*(d-1)/2 + c
                       phasecd = .true.
                    endif

                    if(pair1 >= pair2)then
                      indx = pair1*(pair1-1)/2 + pair2

                    else
                      indx = pair2*(pair2-1)/2 + pair1
                    endif
                     

C------------------------ DECOUPLE MATRIX ELEMENT
                    vtmp = 0.0
                    m = (mpa + mnb)/2
                    jmin = abs(m)
                    jmin = max(jmin,vtbme(indx)%jmin)
                    Jmax = min( jpa + jnb, jpc+jnd)/2
                    jmax = min(jmax,vtbme(indx)%jmax)

                    if(jmin > jmax)cycle


                    do J = Jmin,Jmax
                      fact0 = 1
                      fact1 = 1
                      if(phaseab)then
                        fact0 = (-1)**((jpa+jnb)/2+J)
                        fact1 = - fact0
                      endif
                      if(phasecd)then
                        fact0 = fact0*(-1)**((jpc+jnd)/2+J)
                        fact1 = -fact1*(-1)**((jpc+jnd)/2+J)
                      endif

                      vtmp = vtmp+ 0.5*(fact0*vtbme(indx)%v(j,0) +
     &                                  fact1*vtbme(indx)%v(j,1))
     &     *zeta(a,b)*zeta(c,d)
     &     *cleb(jpa,mpa,jnb,mnb,2*j,2*m)*cleb(jpc,mpc,jnd,mnd,2*j,2*m)

                    enddo  ! loop over J
                    if(abs(vtmp) < 0.00001)cycle
                    nct = nct + 1
                    if(.not. countflag)then
                       hmatPN(nct) = vtmp
                       hmatorbPN(1,nct) = pa
                       hmatorbPN(2,nct) = nb
                       hmatorbPN(3,nct) = pc
                       hmatorbPN(4,nct) = nd

                    endif

                  enddo ! loop over nd
             enddo  ! loop over nb

          enddo  ! loop over pc
       enddo  ! loop over pa
       if(countflag)then
           nmatPN = nct
           print*,' There are ',nmatPN,' PN matrix elements '
           if(allocated(hmatPN))then
                  deallocate(hmatPN)
           end if
           allocate(hmatPN(nct))
           if(allocated(hmatorbPN))then
                  deallocate(hmatorbPN)
           end if
           allocate(hmatorbPN(4,nct))
       endif

       return
       end subroutine countcreatetbmePN

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       real function zeta(i,j)

      implicit none

      integer i,j

      zeta = 1.0
      if(i ==j)zeta = sqrt(2.)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine undoSPE
C
C  at this time assume spes diagonal however can generalize
C  allow assume neutron, proton part the same

      use sporbit
      use spstate
      use interaction

      integer asp,a
      integer it
      if(allocated(speunX)) then
            deallocate(speunX)
      end if

      allocate(speunX(nsps,nsps,2))

      speunX = 0.0
      do it = 1,2
        do asp = 1,nsps
           a = spsqn(asp)%orb
           speunX(asp,asp,it) = spe(a)
        enddo

      enddo
      return
      end subroutine undoSPE

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine undoSPEshifted
C
C  at this time assume spes diagonal however can generalize
C  allow assume neutron, proton part the same

      use sporbit
      use spstate
      use interaction

      integer asp,a
      integer it

      allocate(speunXshift(nsps,nsps,2))

      speunXshift = 0.0
      do it = 1,2
        do asp = 1,nsps
           a = spsqn(asp)%orb
           speunXshift(asp,asp,it) = spe(a)+speshift(a)
        enddo

      enddo
      return
      end subroutine undoSPEshifted


