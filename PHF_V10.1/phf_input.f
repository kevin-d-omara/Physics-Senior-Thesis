CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine get_sp_info
C
C  reads in s.p. space information;
C  fills in q#s into orbqn and spsqn
C
C  reads in either .spo or.sps files
C  it first looks for .spo file; if that fails, 
C  automatically looks for .sps file of same name
C
C----------------------------------------------

      use sporbit
      use spstate
      implicit none

C------ FILE CONTROL ---------------------------

      character*25 filename
      character*1 achar
       
      character*70 title
      integer ilast


C------------ TEMP -------------

      real :: xn,xl,xj,xw    ! orbital q#s
C-----------------------------------------------
C      dummy counters
C--------------------------------------------------
      integer i,j,k,m
      logical success
      integer isp
C---------------BEGIN -------------------------

C-------------- OPEN A FILE ----------------------------------
      success = .false.

      do while(.not.success)
          print*,' Enter file with s.p. orbit information (.spo/.sps)'
          read(5,'(a)')filename
          ilast = index(filename,' ')-1
C........... ATTEMPT TO OPEN .spo FILE.................
          open(unit=1,file=filename(1:ilast)//'.spo',status='old', 
     &     err=101)
          success = .true.
          cycle
101       continue
C...........ATTEMPT TO OPEN .sps FILE..................
          open(unit=1,file=filename(1:ilast)//'.sps',status='old',
     &      err=102)
          success = .true.
          cycle
102       continue
          print*,filename(1:ilast),'.spo/.sps file does not exist '

      enddo
C-------------- READ PAST TITLE CARDS ---------------------------

      success = .false.
      do while(.not.success)
        read(1,'(a)')achar
        backspace(1)
        if(achar /= '#' )then
           success = .true.
        else
           read(1,'(a)')title
           write(6,*)title
        endif
      enddo

C-------------- READ PAST POSSIBLE LABEL OF ISO/PN

      read(1,'(a)')achar
      if(achar == 'p' .or. achar=='P')then
          print*,' .sps file in pn formalism, cannot handle '
          stop
      elseif(achar /= 'i' .and. achar/='I')then
          backspace(1)
      endif

C............ READ # OF ORBITS----------

      read(1,*)numorb

C----------------ALLOCATE MEMORY ------------
      if(numruns.ne.0) then
            deallocate(orbqn,spsqn)
      end if

      allocate(orbqn(numorb))

C---------------READ IN---------------------

      isp = 1
      do i = 1,numorb
        read(1,*,end=2001)xn,xl,xj,xw
        orbqn(i)%nr = int(xn)
        orbqn(i)%l = int(xl)
        orbqn(i)%j = int(2*xj)
        orbqn(i)%w = int(xw)
        orbqn(i)%spstart = isp
        isp = isp + orbqn(i)%j+1
      enddo

      close(unit=1)

C------------ SET UP S.P. STATE INFO -------------

      nsps = 0
      do i = 1,numorb
          nsps = nsps+orbqn(i)%j+1
      enddo

      allocate(spsqn(nsps))
      k = 0
      do i = 1,numorb
        j = orbqn(i)%j
        do m = -j,j,2
           k = k +1
           spsqn(k)%nr = orbqn(i)%nr
           spsqn(k)%l = orbqn(i)%l
           spsqn(k)%par = (-1)**(orbqn(i)%l)
           spsqn(k)%j = orbqn(i)%j
           spsqn(k)%w = orbqn(i)%w
           spsqn(k)%m = m
           spsqn(k)%orb = i
        enddo
      enddo

      return
C------------- ERROR TRAP FOR END OF FILE -----------
2001  continue
      print*,' sudden end of file in ',filename(1:ilast)
      stop

      end subroutine get_sp_info

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine setup4tbmes

C
C   V(a b, c d; J T)
C   assume:  a >= b, c >= d,  a >= c etc.

      use sporbit


      use interaction
      implicit none
     
      integer iorb
      integer icount

      integer npair
      integer pair1,pair2
      integer indx
      integer jmin,jmax

      integer ia,ib,ic,id
      integer na,nb,nc,nd
      integer ja,jb,jc,jd
      integer dstart



C      SOME INTRONS

      norb_allow = numorb
      if (numruns.ne.0) then
            deallocate(orblist,spe,vtbme)
      end if
      allocate(orblist(norb_allow))
      allocate(spe(numorb))

      icount = 0

      do iorb = 1,numorb
           icount = icount+1
           orblist(icount) = iorb
      enddo
C------------ figure out approx dimension of unique matrix elements

      npair = norb_allow*(norb_allow+1)/2  ! # of pairs

      nvdim = npair*(npair+1)/2
      allocate(vtbme(nvdim))

C---------------FIND MIN,MAX J

      do ia = 1,norb_allow
        na = orblist(ia)
        ja = orbqn(na)%j
        do ib = 1,ia
           nb = orblist(ib)
           jb = orbqn(nb)%j
           pair1 = ia*(ia-1)/2 + ib
           do ic = 1,ia
             nc = orblist(ic)
             jc = orbqn(nc)%j

             if(ia ==ic)then
                dstart = ib
             else
                dstart = ic
             endif
             do id = 1,dstart 
                nd = orblist(id)
                jd = orbqn(nd)%j
                pair2 = ic*(ic-1)/2+id
                indx = pair1*(pair1-1)/2+pair2
C                write(18,202)indx,ia,ib,ic,id,pair1,pair2
202         format(i3,3x,4i2,3x,2i2)
                jmax = min((ja+jb)/2, (jc+jd)/2)
                jmin = max(abs(ja-jb)/2, abs(jc-jd)/2)
                vtbme(indx)%jmax = jmax
                vtbme(indx)%jmin = jmin
                
                if(jmin <= jmax)then
                  allocate(vtbme(indx)%v(jmin:jmax,0:1))
                  vtbme(indx)%v = 0.0
                endif

C------------------ STORE ORBITAL INDICES 

                vtbme(indx)%orb(1) = ia
                vtbme(indx)%orb(2) = ib
                vtbme(indx)%orb(3) = ic
                vtbme(indx)%orb(4) = id


             enddo  ! loop over id

           enddo  ! loop over ic

        enddo  ! loop over ia

      enddo  ! loop over ib

      end subroutine setup4tbmes

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine readvtbme
C
C
C

      use sporbit

      use interaction
      implicit none


C------ FILE CONTROL ---------------------------

      character*25 filename
      character*1 achar
       
      character*70 title
      integer ilast

C------------------------------

      integer ia,ib,ic,id
      integer na,nb,nc,nd
      integer pair1,pair2,indx
      integer j,t
      real V
      real, allocatable :: spetmp(:)
      real spscale,a,b,x,vscale  ! for scaling interactions
      integer phase
      integer nme

      integer dw
C-----------------------------------------------
C      dummy counters
C--------------------------------------------------
      integer i,m,mmax
    
      integer L
      logical success
      logical smint   !  a successor to .int 
      logical finished 
C---------------BEGIN -------------------------

      spe = 0.
      finished = .false.
C      print*,maxwtot-minwtot

      do while(.not.finished)

C-------------- OPEN A FILE ----------------------------------
      success = .false.

      do while(.not.success)
          print*,' Enter interaction file name (.smint/.int)'
          print*,' (Enter END to stop ) '
          read(5,'(a)')filename
          if(filename=='END' .or. filename=='end')then
             finished = .true.
             return
          endif
          ilast = index(filename,' ')-1
C........... ATTEMPT TO OPEN .smint FILE.................
          open(unit=1,file=filename(1:ilast)//'.smint',status='old', 
     &     err=101)
          success = .true.
          smint = .true. 
          cycle
101       continue
C...........ATTEMPT TO OPEN .int FILE..................
          open(unit=1,file=filename(1:ilast)//'.int',status='old',
     &      err=102)
          success = .true.
          smint = .false.
          cycle
102       continue
          print*,filename(1:ilast),'.smint/.int file does not exist '

      enddo  ! loop over success

      
C-------------- READ PAST TITLE CARDS ---------------------------

      success = .false.
      do while(.not.success)
        read(1,'(a)')achar
        backspace(1)
        if(achar /= '#' )then
           success = .true.
        else
           read(1,'(a)')title
           write(6,*)title
        endif
      enddo

C--------------- IF THIS IS A .SMINT FILE, COMPARE S.P. STATES 
C       not yet implemented

C-------------- ENTER SCALING -----------------------------

      print*,' Enter scaling for spes, A,B,X ( (A/B)^X ) for TBMEs '
      print*,' (If B or X = 0, then scale by A ) '

      read*,spscale,a,b,x

      if( b == 0. .or. x == 0.0)then
            vscale = a
      else
            vscale = (a/b)**x
      endif
C-------------- READ IN SPEs --------------
      allocate(spetmp(numorb))
      
      read(1,*)nme,(spetmp(i),i= 1,min(10,numorb))

      if(numorb > 10)then
         do m = 10,numorb,10
           mmax = min(10+m,numorb)
           read(1,*)(spetmp(i),i=1+m,mmax)
         enddo
      endif
      spe(:) = spe(:) + spscale*spetmp(:)
          
      deallocate(spetmp)

C-------------- READ IN TBMEs ----------------

      do i = 1,nme

         read(1,*)ia,ib,ic,id,j,t,v
C---------- check for ceiling ----------------

C----------- CHECK PARITY -------------------

         L = orbqn(ia)%l+ orbqn(ib)%l+orbqn(ic)%l+orbqn(id)%l

         if( (-1)**(L) ==-1)then
            print*,' error in parity '
            print*,ia,ib,ic,id,j,t,v
            stop
         endif

C----------- PUT INTO CORRECT ORDER; PICK UP PHASES -------------
C           "CORRECT" ORDER: a >= b, c >= d
        phase = 1
        if(ia < ib)then
           na = ib
           ib = ia
           ia = na
           phase = (-1)**( J+T+(orbqn(ia)%j+orbqn(ib)%j)/2)  ! check
         endif

        if(ic < id)then
           nc = id
           id = ic
           ic = nc
           phase = phase*(-1)**( J+T+(orbqn(ic)%j+orbqn(id)%j)/2)  ! check
         endif

         if(ia < ic .or. (ia==ic .and. ib < id))then
            na = ic
            nb = id
            ic = ia
            id = ib
            ia = na
            ib = nb
         endif

C---------- CONVERT -------------------------
         na = orblist(ia)
         nb = orblist(ib)
         nc = orblist(ic)
         nd = orblist(id)
         pair1 = ia*(ia-1)/2 + ib
         pair2 = ic*(ic-1)/2 + id
         indx = pair1*(pair1-1)/2+pair2
C--------------- ERROR TRAP ---------------

         if(j > vtbme(indx)%jmax .or. j < vtbme(indx)%jmin)then
            print*,' error in Js ',pair1,pair2,indx
            print*,ia,ib,ic,id,J,T,v
            print*,orbqn(ia)%j,orbqn(ib)%j,orbqn(ic)%j,orbqn(id)%j
            print*,vtbme(indx)%jmax,vtbme(indx)%jmin
            
            stop
         endif

         vtbme(indx)%v(j,t)=vtbme(indx)%v(j,t)+v*vscale*phase

      enddo
      close(1)
C---------------------------------------------
      enddo  ! loop over finished
      end subroutine readvtbme

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine readinsd(nsps,n,ifile,psi,errflag)
C
C  subroutine to read in SD from file written by SHERPA
C  NB: probably will have other routines as well
C
C  INPUT:
C  nsps   = # of single-particle states
C  n      = # of particles
C  ifile  = logical number of the file
C  
C  OUTPUT
C   psi(i,j) = slater determinant (real)
C   errflag  = true if there was a problem in reading
C


      implicit none
      integer nsps,n
      real psi(nsps,n) ! slater determinant

      integer i,j
      integer ifile
      logical errflag

      errflag=.false.

      do i = 1,n
C        read(ifile,*,err=103,end=103)(psi(j,i),j=1,nsps)
        read(ifile,err=103,end=103)(psi(j,i),j=1,nsps)
      enddo
      return
  103 continue
      errflag=.true.
      return

      end   
      	

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine GetSherpaSD(pSD,nSD)
C
C  master routine for retrieving SD written out by SHERPA
C  must also convert them from real*4 to complex*8
C
      use system_parameters
      use spstate
      implicit none

      integer ifile
      complex(kind = 8) ::pSD(nsps,numprot),nSD(nsps,numneut)
      real(kind=4),allocatable :: sdtmp(:,:)

      ifile = 73

      call OpenSherpaSD(ifile,nsps,numprot,numneut)

      if(allocated(sdtmp))deallocate(sdtmp)
      if(numprot > 0)then
        allocate(sdtmp(nsps,numprot))

        call ReadSherpa(ifile,nsps,numprot,sdtmp)
        call ConvertSD(nsps,numprot,sdtmp,pSD)

        deallocate(sdtmp)
      endif
      if(numneut > 0)then
        allocate(sdtmp(nsps,numneut))
        call ReadSherpa(ifile,nsps,numneut,sdtmp)
        call ConvertSD(nsps,numneut,sdtmp,nSD)

        deallocate(sdtmp)
      endif
      close(ifile)

      return
      end subroutine GetSherpaSD

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

      subroutine OpenSherpaSD(ifile,nsps,Z,N)

      use sporbit
      implicit none

C.......... 
      integer N,Z
      integer nsps
      integer ifile
C..............FILE HANDLING..........................


      character ychar*1
      character filename*25  		! 
      integer ilast
      integer tempfile			! location of temporary file
      data tempfile/99/
      logical errflag
      character title*60

      logical openflag,failflag

C-------------- DUMMIES 
      integer i,ii,j
      integer norb
      integer zz,nn

C-------------- OPEN FILE ---------------

1     continue

      openflag = .false.
      do while(.not.openflag)

         write(6,*)' Enter input filename (.sd)'
        read(5,'(a)')filename
        ilast=index(filename,' ')
        if(ilast.ne.0)then
          ilast=ilast-1
        else
          ilast=15
        endif

        open(unit=ifile,file=filename(1:ilast)//'.sd',status='old',
     & err=33, form='unformatted')
        openflag = .true.
33      continue

        if(.not.openflag)then
           write(6,*)' That file does not exist; ',
     &    'do you wish to try another file (y/n)?'
           read(5,'(a)')ychar
           if(ychar == 'N' .or. ychar == 'n')stop
        endif

      enddo  ! while on openflag

C..............READ IN HEADER INFO...............
C              CHECK THAT MATCHES SYSTEMS
C

      read(ifile)norb
      failflag = .false.
      if(norb /= numorb)then
	write(6,*)' # of orbits mismatch ',norb,numorb
        failflag = .true.
      endif
      do i = 1,numorb
	read(ifile)ii,norb,j
	if(norb /= orbqn(i)%nr .or. j /= orbqn(i)%j)then
	  write(6,*)' mismatch n,j:',norb,orbqn(i)%nr,j,orbqn(i)%j
	  failflag = .true.
        endif
      enddo

      if(failflag)then

        write(6,*)' The single particle space does not match ',
     &' with that previously chosen.'
        write(6,*)' Exit (x) or choose another file (c)?'
        read(5,'(a)')ychar
        if(ychar.eq.'x' .or. ychar.eq.'X')stop
        close(ifile)
        goto 1
      endif

C...............CHECK IF N,Z match........

      read(ifile)zz,nn
      if(n.ne.nn .or. z.ne.zz)then
	write(6,*)' mismatch Z,N. Old: ',z,n,
     & ', new: ',zz,nn
        write(6,*)' Exit (x) or choose another file (c)?'
        read(5,'(a)')ychar
        if(ychar.eq.'x' .or. ychar.eq.'X')stop
        close(ifile)
        goto 1
      endif

      read(ifile)title

C      ilast=index(title,' ')
C      if(ilast>1)then
C        ilast=ilast-1
        write(6,'(60a)')'Title card: "',title,'"'
C      else
C       print*, 'no title card'
C      endif

      return
      end subroutine OpenSherpaSD

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine ReadSherpa(ifile,nsps,np,sd)

C....... NOTE THIS ASSUMES real*4 
      implicit none

      integer nsps,np
      real(kind=4),intent(OUT) :: sd(nsps,np) ! slater determinant

      integer i,j
      integer ifile
      logical errflag

      errflag=.false.

      do i = 1,np
        read(ifile,err=103,end=103)(sd(j,i),j=1,nsps)

      enddo
      return
  103 continue
      errflag=.true.
      return

      return
      end subroutine ReadSherpa

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      

      subroutine ConvertSD(nsps,np,sd,zsd)

      implicit none
      integer nsps,np

      real(kind=4) :: sd(nsps,np)
      complex(kind = 8) :: zsd(nsps,np)

      integer i,j

      do i = 1,nsps
         do j = 1,np
            zsd(i,j) = dcmplx( sd(i,j),0.0)
         enddo     ! loop over j
      enddo  ! loop over i

      return
      end subroutine ConvertSD

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC      
      subroutine read_sd_txt(psdtrial,nsdtrial)
C=====================================================================
C Reads a (time-reversal) Slater determinant
C in text format
C=====================================================================

      use system_parameters
      use spstate

      implicit none


      integer  :: i,j,k,n,indx,ilast  
      complex(kind=8) :: psdtrial(nsps,numprot),nsdtrial(nsps,numneut)
      real vec(nsps)

      character :: ychar
      character*25 :: filename

  311 continue
      write(6,*)' Enter input filename (.tsd)'
      read(5,'(a)')filename
      ilast=index(filename,' ')
      if(ilast.ne.0)then
          ilast=ilast-1
      else
          ilast=15
      endif

      open(unit=10,file=filename(1:ilast)//'.tsd',status='old',err=33,
     & form='formatted')
      goto 44
  33  continue
      write(6,*)' That file does not exist; ',
     &   'do you wish to try another file (y/n)?'
      read(5,'(a)')ychar
      if(ychar.eq.'n' .or.ychar.eq.'N')then
         return
      else
        goto 311
      endif
  44  continue

      psdtrial(:,:)=dcmplx(0.d0,0.d0)
      if(numneut/=0)nsdtrial(:,:)=dcmplx(0.d0,0.d0)

      do n=1,numprot
          read(10,*)(vec(i),i=1,nsps)
          do i=1,nsps
             psdtrial(i,n)=dcmplx(dble(vec(i)),0.d0)
          enddo
       enddo
       do n=1,numneut
          read(10,*)(vec(i),i=1,nsps)
          do i=1,nsps
             nsdtrial(i,n)=dcmplx(dble(vec(i)),0.d0)
          enddo

      enddo
  122 continue
      close(10)

      return
      end subroutine read_sd_txt
