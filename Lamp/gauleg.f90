SUBROUTINE gauleg(x1,x2,x,w,n)

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL (KIND = 8), INTENT(IN) :: x1, x2
REAL (KIND = 8), DIMENSION(1:n),INTENT(OUT) :: x, w

INTEGER :: i,j,m
REAL (KIND = 8) :: EPS, p1, p2, p3, pp, xl, xm, z, z1

PARAMETER(EPS = 3.0d-14)

m = (n+1)/2
xm = 0.5d0*(x2+x1)
xl = 0.5d0*(x2-x1)

DO i = 1,m
	z = DCOS(3.1415926535897932d0*(DBLE(i)-0.25d0)/(DBLE(n)+0.5d0))
1	CONTINUE
		p1  = 1.0d0
		p2 = 0.0d0
		DO j = 1,n
			p3 = p2
			p2 = p1
			p1 = ((2.0d0*DBLE(j) - 1.0d0)*z*p2 - (DBLE(j)-1.0d0)*p3)/DBLE(j)
		END DO
		pp = n*(z*p1-p2)/(z*z - 1.0d0)
		z1 = z
		z = z1 - p1/pp
	IF (DABS(z-z1) > EPS) GOTO 1
	x(i) = xm - xl*z
	x(n+1 - i) = xm + xl*z
	w(i) = 2.0d0*xl/((1.0d0-z*z)*pp*pp)
	w(n+1-i) = w(i)
END DO
RETURN

END SUBROUTINE gauleg