REAL (KIND = 8) FUNCTION gammln(xx)			!j is a double, use DBLE(n) to convert when sending variable
							!g_fact is a double
IMPLICIT NONE

INTEGER :: j
REAL (KIND = 8) :: xx,ser,stp,tmp,x,y,cof(6)

SAVE cof, stp

DATA cof, stp/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,-1.231739572450155d0,0.1208650973866179d-2,&
			-0.5395239384953d-5,2.5066282746310005d0/

x=xx
y=x

tmp = x+5.5d0
tmp = (x+0.5d0)*DLOG(tmp) - tmp
ser = 1.000000000190015d0

DO j = 1,6
	y = y+1.0d0
	ser = ser+cof(j)/y
END DO

gammln = tmp +DLOG(stp*ser/x)

END FUNCTION gammln
      

SUBROUTINE gaujac(x,w,n,alf,bet)
      
IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL (KIND = 8), INTENT(IN) :: alf, bet
REAL (KIND = 8), INTENT(INOUT) :: x(n), w(n)

INTEGER ::MAXIT
REAL (KIND = 8) :: EPS

PARAMETER (EPS=3.0D-14,MAXIT=100)
INTEGER i,its,j
REAL (KIND = 8) :: alfbet,an,bn,r1,r2,r3,gammln
REAL (KIND = 8) :: a,b,c,p1,p2,p3,pp,temp,z,z1
      do 13 i=1,n
        if(i == 1)then
          an=alf/n
          bn=bet/n
          r1=(1.0d0+alf)*(2.78d0/(4.d0+DBLE(n*n))+0.768d0*an/DBLE(n))
          r2=1.0d0+1.48d0*an+0.96d0*bn+0.452d0*an*an+0.83d0*an*bn
          z=1.0d0-r1/r2
        else if(i == 2)then
          r1=(4.1d0+alf)/((1.+alf)*(1.0d0+0.156d0*alf))
          r2=1.+.06*(DBLE(n)-8.d0)*(1.0d0+0.12d0*alf)/DBLE(n)
          r3=1.0d0+0.012d0*bet*(1.0d0+0.25d0*dabs(alf))/DBLE(n)
          z=z-(1.0d0-z)*r1*r2*r3
        else if(i == 3)then
          r1=(1.67d0+.28d0*alf)/(1.0d0+0.37d0*alf)
          r2=1.0d0+0.22d0*(DBLE(n)-8.0d0)/DBLE(n)
          r3=1.0d0+8.0d0*bet/((6.28d0+bet)*DBLE(n*n))
          z=z-(x(1)-z)*r1*r2*r3
        else if(i==n-1)then
          r1=(1.0d0+0.235d0*bet)/(0.766d0+0.119d0*bet)
          r2=1.0d0/(1.0d0+0.639d0*(DBLE(n)-4.0d0)/(1.0d0+0.71d0*(DBLE(n)-4.0d0)))
          r3=1.0d0/(1.0d0+20.0d0*alf/((7.5d0+alf)*DBLE(n*n)))
          z=z+(z-x(n-3))*r1*r2*r3
        else if(i==n)then
          r1=(1.0d0+0.37d0*bet)/(1.67d0+0.28d0*bet)
          r2=1.0d0/(1.0d0+0.22d0*(DBLE(n)-8.0d0)/DBLE(n))
          r3=1.0d0/(1.0d0+8.0d0*alf/((6.28d0+alf)*DBLE(n*n)))
          z=z+(z-x(n-2))*r1*r2*r3
        else
          z=3.0d0*x(i-1)-3.0d0*x(i-2)+x(i-3)
        endif
        alfbet=alf+bet
        do 12 its=1,MAXIT
          temp=2.0d0+alfbet
          p1=(alf-bet+temp*z)/2.0d0
          p2=1.0d0
          do 11 j=2,n
            p3=p2
            p2=p1
            temp=2*j+alfbet
            a=2*j*(j+alfbet)*(temp-2.0d0)
            b=(temp-1.0d0)*(alf*alf-bet*bet+temp*(temp-2.0d0)*z)
            c=2.0d0*(j-1+alf)*(j-1+bet)*temp
            p1=(b*p2-c*p3)/a
11        continue
          pp=(n*(alf-bet-temp*z)*p1+2.0d0*(DBLE(n)+alf)*(DBLE(n)+bet)*p2)/(temp*(1.0d0-z*z))
          z1=z
          z=z1-p1/pp
          if(dabs(z-z1)<=EPS)goto 1
12      continue
        write(*,*) 'too many iterations in gaujac'
        read(*,*)  !pause until user hits enter
1       x(i)=z
        w(i)=dexp(gammln(alf+DBLE(n))+gammln(bet+DBLE(n))-gammln(DBLE(n)+1.0d0)-gammln(DBLE(n)+alfbet+1.))&
        		*temp*2.0d0**alfbet/(pp*p2)
	!write(*,*) x(i),w(i)
13    continue
      return
      END
