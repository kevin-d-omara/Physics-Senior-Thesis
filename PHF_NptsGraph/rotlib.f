CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  package MONTY_ROTLIB
C
C  routines for rotation
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine testrotate
C
C  master subroutine to test rotations and others
C
      use system_parameters

      use spstate
      implicit none


      complex(kind = 8),allocatable :: psdr(:,:),nsdr(:,:)
      complex(kind = 8),allocatable :: psdl(:,:),nsdl(:,:)

      complex(kind = 8),allocatable :: rhopij(:,:),rhonij(:,:)
      complex(kind = 8) :: ovlpp,ovlpn

      complex(kind = 8) :: vme

      real :: a  ! angle of rotation
      complex(kind = 8),allocatable :: drotZ(:,:)
      integer ia


      allocate(psdr(nsps,numprot),nsdr(nsps,numneut))
      allocate(psdl(nsps,numprot),nsdl(nsps,numneut))

      call GetSherpaSD(psdr,nsdr)

      allocate(rhopij(nsps,nsps),rhonij(nsps,nsps))
      call makerhoij(1,numprot,psdr,psdr,ovlpp,rhopij)
      call makerhoij(2,numneut,nsdr,nsdr,ovlpn,rhonij)

      print*,' Overlaps ',ovlpp,ovlpn
C      write(17,*)rhopij

      call TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
      print*,' < H > = ',vme

C      print*,' '
C      print*,' Enter angle of rotation in degrees '
C      read*,a
C      a = a*3.1415926/180.

      allocate(drotZ(nsps,nsps))

      do ia = 1,360
      a = ia*3.1415926/180.

      call makeZrot(nsps,a,drotZ)
      print*,' rotation made '
      call rotate(nsps,numprot,drotZ,psdr,psdl)
      call rotate(nsps,numneut,drotZ,nsdr,nsdl)
C      print*,' rotated!'
      call makerhoij(1,numprot,psdl,psdr,ovlpp,rhopij)
      call makerhoij(2,numneut,nsdl,nsdr,ovlpn,rhonij)

      print*,' Overlaps ',ovlpp,ovlpn
C      write(17,*)rhopij

      call TBMEmaster(nsps,rhopij,rhonij,ovlpp,ovlpn,vme)
      print*,' < H > = ',vme
      write(17,*)ia,vme/ovlpp/ovlpn
      enddo
      return

      end subroutine testrotate


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine makeZrot(ndim,a,drotZ)
C
C  creates matrix to rotates about z-axis
C  exp( i a Jz )
C

      use spstate
      implicit none
      integer ndim
      complex(kind = 8) :: drotZ(ndim,ndim)

      real a
      real(kind = 8) :: xm
      complex(kind = 8) :: i,zzero

      integer m,n

      zzero = (0.d0,0.d0)
      i     = (0.d0,1.d0)

C      drotZ(:,:) = zzero

      if(nsps > ndim)then
         print*,' mismatch in nsps ndim ',nsps,ndim
         stop
      endif
      do n = 1,nsps
        do m = 1,nsps
          drotZ(m,n) = zzero
        enddo
        xm = dfloat( spsqn(n)%m )*0.5d0
        drotZ(n,n) = cmplx(dcos(a*xm),dsin(a*xm))
C        print*,n,drotZ(n,n)
      enddo

      return
      end subroutine makeZrot

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine rotate(n,npart,drot,psiin,psiout)
      implicit none

      integer n
      integer npart
      complex(kind = 8) :: drot(n,n)
      complex(kind = 8) :: psiin(n,npart),psiout(n,npart)

      integer i,j,k
      complex(kind = 8) :: ztmp

C      psiout = (0.d0,0.d0)

      do j = 1,npart
         do i = 1,n
             ztmp = (0.d0,0.d0)
             do k = 1,n
                ztmp = ztmp + drot(i,k)*psiin(k,j)
             enddo
             psiout(i,j) = ztmp

         enddo  ! loop over j
      enddo ! loop over i

      return
      end subroutine rotate

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C
C    routines for computing wigner d-function
C    originals by WEO, modified by CWJ  3/98  LSU
C
      subroutine makeYrot(ndim,a,drotY)
C
C  routine to create a rotation matrix
C	ASSUME PROTONS = NEUTRONS spaces
C
C  INPUT:
C	nsps = # of single-particle m-states
C	spsqn    quantum numbers of m-states
C	angle = angle of rotation
C
C  OUTPUT:
C	rotate = matrix of rotation
C
C  FUNCTIONS CALLED
C	wigner_d
C
      use spstate
      implicit none
      integer ndim
      complex(kind = 8) :: drotY(ndim,ndim)

      real a
      
      real wigner_d	! wigner d-function
      
      real xj,xm,xmp
      
      integer k,l	! dummy indices
      if(nsps > ndim)then
         print*,' Error in rotation Y '
         stop
      endif
      do k = 1,nsps
      	xj = float(spsqn(k)%j)*0.5
      	xmp = float(spsqn(k)%m)*0.5
      	do l = 1,nsps
      	    if( spsqn(k)%orb == spsqn(l)%orb )then
      	    	xm = float(spsqn(l)%m)*0.5
      	    	drotY(k,l) = dcmplx(wigner_d(xj,xmp,xm,a),0.0)
      	    else
      	    	drotY(k,l)= (0.d0,0.d0)
      	    endif
      	enddo
      enddo 

      return 
      end subroutine makeYrot

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      function wigner_d(xjj,xmp,xm,theta)

C------------computes Wigner (little) d-function d^J_m'm(beta)
C------------cf Edmonds 
C     in double precision

      implicit none

      real(kind = 4) wigner_d
C---------INPUT--------------------------

      real xjj,xmp,xm		! coefficients of d^j_m'm
      real theta		! angle of rotation 

C---------INTERMEDIARIES-----------------
      
      integer i1,i2,i3,i4	! various combination of j-m etc
      real alpha,beta		! alpha = m'-m,  beta = m'+m
      real phase		! from symmetries of d-function 
      				! must have alpha, beta > 0 
      real xmm,xmmp 		! from reordering of m',m 
      real uu			! argument = cos(theta)
      integer n			! order of jacobi polynomial = j-m'
      real xjacobi		! function = value jacobi polynomial
      real xnorm		! overall factor
      
      integer size		! dimension for factorial array 
      parameter(size=40)
      real fac_ar(0:size)
      integer imax

C-------------------------------------------------------------

      if(abs(xmp).gt.xjj)stop 'xmp > xjj in wigner_d '
      if(abs(xm).gt.xjj)stop 'xm > xjj in wigner_d '
      
c------------------------    put M values into order for
c                            Jacobi polynomials

      call jac_param(xjj,xmp,xm,xmm,xmmp,n,alpha,beta,phase)

      i1=nint(xjj+xmm)
      i2=nint(xjj-xmm)
      i3=nint(xjj+xmmp)
      i4=nint(xjj-xmmp)
     
      imax = max(i1,i2)
      imax = max(imax,i3)
      imax = max(imax,i4)
      call lnfact(imax,size,fac_ar)

      xnorm=exp((fac_ar(i1)+fac_ar(i2)-fac_ar(i3)-fac_ar(i4))/2.)

      uu=cos(theta)
      i1=nint(beta)
      i2=nint(alpha)

      wigner_d=phase*xnorm*cos(theta/2.)**i1*sin(theta/2.)**i2*
     &        xjacobi(n,alpha,beta,uu)

      return
      end function wigner_d

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine lnfact(imax,size,fac_ar)
C
C  returns array fac_ar of ln factorials
C      
      integer imax		! max factorial needed
      integer size
      real fac_ar(0:size)	! fac_ar(i) = ln(i!)
      integer i
      
      if(imax.gt.size)stop ' error in lnfact '
      fac_ar(0) = 0.0
      fac_ar(1) = 0.0
      do i = 2,imax
        fac_ar(i)=fac_ar(i-1)+log(float(i))
      enddo
      return
      end
      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine jac_param(xjj,xm,xmp,xmm,xmmp,
     &                     n,alpha,beta,phase)

C------------- returns parameters and, if necessary, a phase for 
C              input into jacobi polynomials

      implicit none 

C------------INPUT------------------------------------------

      real xjj,xm,xmp		! from wigner d^j_m'm
      
C-----------OUTPUT------------------------------------------

      real xmm,xmmp		! reordering/symmetry of m',m
      real alpha,beta		! alpha, beta >= 0
      				! alpha = m'-m,  beta = m'+m
      real phase		! from symmetries/reordering
      integer n			! order of polynomial = j - m'				

C-----------INTERMEDIARIES

      integer iph

C----------------------------------------------------------      
      n=0
      alpha=0.0
      beta=0.0
      phase=1.0

      xmm = xm
      xmmp = xmp
      if(xjj.eq.0.00)return
      if(xm.ge.0.00)then
         if(xmp.ge.0.00)then
             if(xm.ge.xmp)then
                xmm=xm
                xmmp=xmp
                phase=1.0
             else
                xmm=xmp
                xmmp=xm
                iph=nint(xm-xmp)
                phase=(-1.0)**iph
             end if
         elseif(xmp.lt.0.0)then
             if(xm.ge.abs(xmp))then
                xmm=xm
                xmmp=xmp
                phase=1.00
             else
                xmm=-xmp
                xmmp=-xm
                phase=1.0
             end if
         end if
      elseif(xm.lt.0.0)then
         if(xmp.ge.0.0)then
             if(abs(xm).ge.xmp)then
                xmm=-xm
                xmmp=-xmp
                iph=nint(xm-xmp)
                phase=(-1.0)**iph
             else
                xmm=xmp
                xmmp=xm
                iph=nint(xm-xmp)
                phase=(-1.00)**iph
             end if
         elseif(xmp.lt.0.0)then
             if(abs(xm).ge.abs(xmp))then
                xmm=-xm
                xmmp=-xmp
                iph=nint(xm-xmp)
                phase=(-1.0)**iph
             else
                xmm=-xmp
                xmmp=-xm
                phase=1.0
             end if
         end if
      end if

      alpha=xmm-xmmp
      beta=xmm+xmmp
      if(alpha.lt.0.0)then
          write(6,*)xm,xmp,xmm,xmmp
          stop 'error in wigner_d alpha < 0'
      end if
      if(beta.lt.0.0)then
          write(6,*)xm,xmp,xmm,xmmp
         stop 'error in wigner_d alpha < 0'
      end if
      n=nint(xjj-xmm)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      function xjacobi(n,alpha,beta,x) 
C
C  computes jacobi polynomial P^(alpha, beta)_n(x)
C

      implicit none
      real(kind =4) :: xjacobi
c******************   use recurrence relation to compute jacobi polys
      real xjac(0:60)
      integer n		! order of polynomial
      real x		! argument of polynomial
      real alpha,beta

C----------- INTERMEDIATES -----------------------      

      integer i
      real xn
      real one_nab,two_nab
      real a1,a2,a3,a4

C-------------------------------------------------
      
      xjac(0)=1.0
      xjac(1)=0.50*(2.00*(alpha+1.00)+(alpha+beta+2.00)*
     &        (x-1.00))
      do i=2,n
         xn=float(i-1)
         one_nab=xn+alpha+beta
         two_nab=2.0*xn+alpha+beta
         a1=2.0*(xn+1.0)*(one_nab+1.0)*two_nab
         a2=(two_nab+1.0)*(alpha**2-beta**2)
         a3=two_nab*(two_nab+1.0)*(two_nab+2.0)
         a4=2.0*(xn+alpha)*(xn+beta)*(two_nab+2.0)
         xjac(i)=((a2+a3*x)*xjac(i-1)-a4*xjac(i-2))/a1
      end do 
      xjacobi=xjac(n)
      return 
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine make_iJy(zijy)

      use spstate

      implicit none
      complex(kind = 8) :: zijy(nsps,nsps)
      integer ia,ib
      integer ja,ma,la,na
      integer jb,mb,lb,nb

      zijy(:,:) = (0.d0,0.d0)

C---------------------------- 
      do ia = 1,nsps

        ja = spsqn(ia)%j
        ma = spsqn(ia)%m
        la = 2*spsqn(ia)%l
        na = spsqn(ia)%nr
        do ib = 1,nsps
           jb = spsqn(ib)%j
           mb = spsqn(ib)%m
           lb = 2*spsqn(ib)%l
           nb = spsqn(ib)%nr
           if( na /= nb .or. la/=lb .or. ja /=jb)cycle
           if( ma == mb+2) then
            zijy(ia,ib)=dcmplx(0.5*sqrt(float((jb-mb)*(jb+mb+2))/4.))
           endif
           if( ma == mb-2)then
           zijy(ia,ib)=dcmplx(-0.5*sqrt(float((jb+mb)*(jb-mb+2))/4.))
           endif
         enddo
      enddo

      return

      end subroutine make_iJy
                                                                  

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC