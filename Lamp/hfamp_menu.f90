!====================================================================

SUBROUTINE J_Single(psd,nsd,pairTest)

USE system_parameters
USE spstate
USE psis

IMPLICIT NONE
interface
 SUBROUTINE Projection_with_Parity(js, np, psd, nsd, nf, pevals,posdef,PairTest,ntrace)  ! INTERFACE
 USE phf_vals
 USE system_parameters
 USE spstate
 USE psis
 USE errortests
 IMPLICIT NONE
 INTEGER (KIND = 8), INTENT(IN) :: np
 REAL (KIND = 8), INTENT(IN) :: js
 COMPLEX (KIND = 8), INTENT(IN) :: psd(numsd,nsps,numprot), nsd(numsd,nsps,numneut)
 LOGICAL, INTENT(IN) :: PairTest
 LOGICAL, INTENT(INOUT) :: posdef(2)
 INTEGER (KIND = 8), INTENT(OUT) :: nf
 REAL (KIND = 8), INTENT(OUT) :: pevals(npair,numsd*np)
 real (kind=8) :: ntrace(2)
 
end subroutine Projection_with_Parity   ! INTERFACE
end interface

COMPLEX (KIND = 8), INTENT(IN) :: psd(nsps,numprot), nsd(nsps,numneut)
LOGICAL, INTENT(IN) :: pairTest

INTEGER :: A
INTEGER :: clock_start, clock_end, clock_rate
REAL :: elapsed_time
INTEGER (KIND = 8) :: np,nf,i
REAL (KIND = 8) :: js,check
REAL (KIND = 8), ALLOCATABLE, DIMENSION(:) :: pevals, jvals
REAL (KIND = 8), ALLOCATABLE, DIMENSION(:,:) ::PairEvals
CHARACTER (LEN = 1) :: choice
CHARACTER (LEN = 12) :: jstring
LOGICAL :: posdef(2) = .true.


1001 FORMAT(1F8.1,G15.8)
1002 FORMAT(1F8.1,G15.8,G15.8)

A = numprot + numneut
IF (MOD(A,2) == 0) THEN
	check = 0.0d0
	jstring = 'integer'
ELSE
	check = 1.0d0
	jstring = 'half-integer'
ENDIF

WRITE(*,*) 'Enter coupled J value (',TRIM(jstring),' value):'
DO
	READ(*,*) js
	IF (js < 0.0d0) THEN
		WRITE(*,*) 'J cannot be negative, re-enter.'
	ELSEIF (MOD(2.0d0*js,2.0d0) == check) THEN
		EXIT
	ELSE
		WRITE(*,*) 'Incorrect J.  A = ', A
		WRITE(*,*) 'Expecting ', TRIM(jstring), ' values for J'
		WRITE(*,*) 'Please re-enter J.'
	END IF
END DO

np = INT(2.0d0*js) + 1

IF (.NOT.pairTest) THEN
	ALLOCATE(pevals(np))
	CALL SYSTEM_CLOCK(COUNT_RATE = clock_rate)
	CALL SYSTEM_CLOCK(COUNT = clock_start)
	CALL JackedTest(js,np,psd,nsd,nf,pevals,posdef)
	CALL SYSTEM_CLOCK(COUNT = clock_end)
ELSE
	ALLOCATE(PairEvals(2,np))
	CALL SYSTEM_CLOCK(COUNT_RATE = clock_rate)
	CALL SYSTEM_CLOCK(COUNT = clock_start)
!	CALL Projection_with_Parity(js,np,psd,nsd,nf,PairEvals,posdef,pairtest)
	CALL SYSTEM_CLOCK(COUNT = clock_end)
END IF

elapsed_time = (REAL(clock_end) - REAL(clock_start))/REAL(clock_rate)

1099 format(2x,'J',4x,2(G15.8,4x))

!WRITE(*,*) ''
IF ((nf == 0).OR.(posdef(1) .EQV. .false.)) THEN
	RETURN
ELSE
	WRITE(*,*) '    J  ',' Energy (MeV)'
	write(*,*) 'States tot: ',nf
	IF (pairTest) THEN
		WRITE(*,1099)
	END IF
	DO i = 1, nf
		IF (.NOT.pairTest) THEN
			WRITE(*,1001) js,pevals(i)
		ELSE
			WRITE(*,1002) js,PairEvals(1,i),PairEvals(2,i)
		END IF
	END DO
ENDIF
WRITE(*,*)
WRITE(*,*) 'Total projection time: ',elapsed_time
WRITE(*,*)

ALLOCATE(jvals(np))

jvals = DBLE(nf)

CALL J_Write(np,nf,jvals,pevals)

END SUBROUTINE J_Single
!====================================================================
!
!  INPUT:
!    psdf = proton Slater determinant(s)
!    nsdf = neutron Slater determinant(s)
!    pairTest = logical for cases with mixed parity
!
!  CALLED BY:
!    --- NO LONGER USED
!
!  SUBROUTINES CALLED:
!   phf_dealloc
!   phf_alloc
!   JackedTestPair
!   J_WritePair
!   J_Write
!
SUBROUTINE J_Range(psdf,nsdf,pairTest)

USE phf_vals
USE phf_allocation
USE system_parameters
USE spstate
USE psis

IMPLICIT NONE

interface
	
	SUBROUTINE J_Write(np,nf,jvals,pevals)   ! INTERFACE
	USE errortests
	IMPLICIT NONE
	INTEGER (KIND = 8), INTENT(IN) :: nf, np
	REAL (KIND = 8), DIMENSION(np), INTENT(IN) :: jvals
	REAL (KIND = 8), DIMENSION(np), INTENT(IN) :: pevals
    end subroutine J_write   ! INTERFACE

    SUBROUTINE Projection_with_Parity(js, np, psd, nsd, tol,npts,nf, pevals,posdef,PairTest,ntrace)  ! INTERFACE
    USE phf_vals
    USE system_parameters
    USE spstate
    USE psis
    USE errortests
    IMPLICIT NONE
    INTEGER (KIND = 8), INTENT(IN) :: np
    REAL (KIND = 8), INTENT(IN) :: js
    COMPLEX (KIND = 8), INTENT(IN) :: psd(numsd,nsps,numprot), nsd(numsd,nsps,numneut)
	real (kind=4) :: tol
	integer(kind=8):: npts
    LOGICAL, INTENT(IN) :: PairTest
    LOGICAL, INTENT(INOUT) :: posdef(2)
    INTEGER (KIND = 8), INTENT(OUT) :: nf
    REAL (KIND = 8), INTENT(OUT) :: pevals(npair,numsd*np)
	real (kind=8) :: ntrace(2)
   end subroutine Projection_with_Parity   ! INTERFACE
   
end interface

COMPLEX (KIND = 8),INTENT(IN) :: psdf(numsd,nsps,numprot),nsdf(numsd,nsps,numneut)
LOGICAL,INTENT(INOUT) :: pairTest
real (kind=8) :: ntrace(2)

real(kind=8),allocatable :: problist(:,:),jlist(:)  ! probabilities for finding a given J
real(kind=8) :: probsum                    ! summed probability
integer :: nlist,ilist
INTEGER :: A
INTEGER (KIND = 8):: i, k, nf, np, npmax, nftotal
INTEGER :: clock_start, clock_end, clock_rate
REAL :: elapsed_time
REAL(KIND = 8) :: jmin,jmax,js,check,hdtot,ndtot,hd,nd
REAL(KIND = 8), ALLOCATABLE, DIMENSION(:) :: pevals,jall,pall
REAL(KIND = 8), ALLOCATABLE, DIMENSION(:,:) :: PairEvals,pallPair
CHARACTER (LEN = 13) :: jstring
LOGICAL :: posdef = .true.,posdefP(2) = .true.

logical :: ask4tol = .true. 
real :: tol
integer(kind=8) :: npts

1001 FORMAT(I6,4X,F12.5,4X,F8.5,4X,F3.1,4x)
!1001 FORMAT(1F8.1,4x,G15.8)
1002 FORMAT(I6,4x,2(F13.5),4x,F3.1)
1003 FORMAT(3x,'State #    E (MeV)       J')
1033 FORMAT(3x,'State #    E+(MeV)     E-(MeV)       J')

1004 FORMAT(5x,'------------------')
1005 FORMAT(18x,'E(MeV) w/ Parity')
1006 FORMAT(6x,'J',11x,'+',15x,'-')
1007 FORMAT(5x,37('-'))

A = numprot + numneut

IF (MOD(A,2) == 0) THEN
	check = 0.0d0
	jstring = 'integer'
ELSE
	check = 1.0d0
	jstring = 'half-integer'
ENDIF

WRITE(*,*) 'Enter J-min, J-max range to project (',TRIM(jstring),' values):'

DO
	READ(*,*) jmin, jmax
	IF (jmin > jmax) THEN
		WRITE(*,*) 'J-min cannot be larger than J-max. Please Re-enter.'
	ELSEIF (jmin < 0.0d0) THEN
		WRITE(*,*) 'J cannot be negative. Please Re-enter.'
	ELSEIF ((MOD(2.0d0*jmin,2.0d0) == check).AND.(MOD(2.0d0*jmax,2.0d0) == check)) THEN
		EXIT
	ELSE
		WRITE(*,*) 'Incorrect j.  A = ', A
		WRITE(*,*) 'Expecting ', TRIM(jstring), ' values for J'
		WRITE(*,*) 'Please re-enter J-min, J-max.'
	END IF
END DO
WRITE(*,*) ''

js = jmin

CALL phf_dealloc
CALL phf_alloc(jmin,jmax,pairtest)

npmax = INT((jmax-jmin+1.0d0)*(jmin+jmax+1.0d0))

ALLOCATE(pall(npmax),jall(npmax),pallPair(2,npmax))
nlist = int(jmax-jmin+1)
allocate(problist(2, nlist) )
allocate(jlist(nlist) )

problist = 0.d0
jlist= 0.d0
pallpair=0.d0

pall = 0.0d0
jall = 0.0d0

nftotal = 0

!........ ADDED IN V 10..........

if(ask4tol)then
	print*,' Enter tolerance for norm (typical = 0.01) '
	read*,tol
	print*,' Enter # of integration pts (typical = 25)'
	read*,npts
	write(*,*)' # of mesh pts = ',npts,' cutoff criterion = ',tol,''
	
else
	tol= 0.01
	npts = 25
end if 
!...............................
write(*,*)' '
	
!ELSE
	IF (.NOT.pairTest) posdefP(2) = .FALSE.
	if(pairtest)then
		write(*,1033)
	else
		write(*,1003)
	end if
!	WRITE(*,1005)
!	WRITE(*,1006)
	WRITE(*,1007)
	CALL SYSTEM_CLOCK(COUNT_RATE = clock_rate)
	CALL SYSTEM_CLOCK(COUNT = clock_start)
	probsum = 0.0	
	ilist = 0
	DO
		np = INT(2.0d0*js) + 1
		ALLOCATE(PairEvals(npair,numsd*np))		
		CALL Projection_with_Parity(js,np,psdf,nsdf,tol,npts,nf,PairEvals,posdefP(:),pairTest,ntrace)		
		ilist=ilist+1

		if(pairTest)then
		    problist(:,ilist)=ntrace(:)/float(numsd)
			probsum = probsum+ntrace(1)+ntrace(2)
		else
		    problist(1,ilist)=ntrace(1)/float(numsd)
			probsum = probsum + ntrace(1)	/float(numsd)
		end if
		jlist(ilist)=js
		IF ((posdefP(1) .EQV. .TRUE.).OR.(posdefP(2) .EQV. .TRUE.)) THEN
			!WRITE(*,*) ''
			!WRITE(*,*) '    J  ',' Energy (MeV)'
			DO i = 1, nf
				nftotal = nftotal + 1
				IF (pairTest) THEN
					WRITE(*,1002) nftotal,PairEvals(1,i),PairEvals(2,i),js
				ELSE
					WRITE(*,1001) nftotal,PairEvals(1,i),js
				END IF
				jall(nftotal) = js
				pallPair(1,nftotal) = PairEvals(1,i)
				IF (pairTest) pallPair(2,nftotal) = PairEvals(2,i)
			END DO
			!WRITE(*,*) ''
		END IF
		
        if(nf<1)then
			write(*,*)' No states for J = ',js
		end if
		DEALLOCATE(PairEvals)
		js = js + 1.0d0
		IF (js > jmax) THEN
			EXIT
		END IF
	END DO
	CALL SYSTEM_CLOCK(COUNT = CLOCK_END)
	elapsed_time = (REAL(clock_end) - REAL(clock_start))/REAL(clock_rate)
	WRITE(*,*)
	WRITE(*,*) "Total projection time: ",elapsed_time
	WRITE(*,*)
	print*,' Sum of norms = ',probsum
!	IF (pairTest) THEN
!		CALL J_WritePair(npmax,nftotal,jall,pallPair)
!	ELSE
!		print*,npmax,nftotal
!		print*,PairEvals(1,:)
!		CALL J_Write(npmax,nftotal,jall,pallPair(1,:))
!	END IF
!END IF

END SUBROUTINE J_Range

!=============================================================
!	
! INPUT:
!   np = # of possible states
!   nf = # of actual states
!   jvals   : array of Js      (dimension np)
!   pevals : array of energies (dimension np)
! CALLED BY: J_Range
!
!	
SUBROUTINE J_Write(np,nf,jvals,pevals)

USE errortests

IMPLICIT NONE

INTEGER (KIND = 8), INTENT(IN) :: nf, np
REAL (KIND = 8), DIMENSION(np), INTENT(IN) :: jvals
REAL (KIND = 8), DIMENSION(np), INTENT(IN) :: pevals

INTEGER :: i,iostatus
CHARACTER (LEN = 1) :: choice
CHARACTER (LEN = 4) :: shell
CHARACTER (LEN = 26) :: name
CHARACTER (LEN = 100) :: filename

1000 FORMAT(I3,3(G15.8))

WRITE(*,*) 'Write output to file? (Y or N)'

DO
	READ(*,*) choice
	IF ((choice == 'y').OR.(choice == 'Y')) THEN
		EXIT
	ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
		RETURN
	ELSE
		WRITE(*,*) 'Y or N please.'
	ENDIF
END DO

IF (path.EQ."undefined") THEN
	WRITE(*,*) "Enter file path (i.e. ~/Desktop/path/)"
	READ (*,'(a)') path
	!path = "'"//TRIM(path)//"'"
ELSE
	WRITE(*,*) "Current path OK? (Y or N): ",path
	DO
		READ(*,*) choice
		IF ((choice == 'y').OR.(choice == 'Y')) THEN
			EXIT
		ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
			WRITE(*,*) "Enter new path: "
			READ(*,'(a)') path
			path = TRIM(path)
			EXIT
		ELSE
			WRITE(*,*) 'Y or N please.'
		ENDIF
	END DO
END IF

DO
	WRITE(*,*) 'Enter file name (without extention -- .dat added): '
	READ(*,*) name
	filename = TRIM(path)//TRIM(name)//'.dat'
	OPEN(UNIT = 10, FILE = filename, STATUS = 'NEW',IOSTAT = iostatus)
	IF (iostatus > 0) THEN
		WRITE(*,*) "File already exists: Overwrite (o), append (a), new file name (n)?"
		READ(*,*) choice
		DO
			IF ((choice == 'o').OR.(choice == 'O')) THEN
				CLOSE(UNIT = 10)
				OPEN(UNIT = 10, FILE = filename, STATUS = 'REPLACE')
				EXIT
			ELSEIF ((choice == 'a').OR.(choice == 'A')) THEN
				CLOSE(UNIT = 10)
				OPEN(UNIT = 10,FILE = filename, STATUS = 'OLD', POSITION = 'APPEND')
				EXIT
			ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
				CLOSE(UNIT = 10)
				EXIT
			ELSE
				WRITE(*,*) 'Incorrect choice.  Please select again.'
			END IF
		END DO
		IF ((choice.NE.'n').AND.(choice.NE.'N')) EXIT
	ELSE
		EXIT
	END IF
END DO

DO i = 1, nf
	WRITE(10,1000) i,pevals(i),jvals(i)
END DO

WRITE(*,*) ''
WRITE(*,*) 'Data written to:',filename

CLOSE(UNIT = 10)

END SUBROUTINE J_Write
!====================================================================
!
!  revised
! 

SUBROUTINE J_WritePair(np,nf,jvals,pevals)

USE errortests

IMPLICIT NONE

INTEGER (KIND = 8), INTENT(IN) :: nf, np
REAL (KIND = 8), DIMENSION(np), INTENT(IN) :: jvals
REAL (KIND = 8), DIMENSION(2,np), INTENT(IN) :: pevals

INTEGER :: i,iostatus
CHARACTER (LEN = 1) :: choice
CHARACTER (LEN = 4) :: shell
CHARACTER (LEN = 26) :: name
CHARACTER (LEN = 100) :: filename!, path

!path = '/home/Jtstaker/Desktop/JohnsonResearch/ProjectedData/'//TRIM(shell)//'_data/'

1000 FORMAT(I3,4(G15.8))

WRITE(*,*) 'Write output to file? (Y or N)'

DO
	READ(*,*) choice
	IF ((choice == 'y').OR.(choice == 'Y')) THEN
		EXIT
	ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
		RETURN
	ELSE
		WRITE(*,*) 'Y or N please.'
	ENDIF
END DO

IF (path.EQ."undefined") THEN
	WRITE(*,*) "Enter file path (i.e. ~/Desktop/path/)"
	READ (*,'(a)') path
	!path = "'"//TRIM(path)//"'"
ELSE
	WRITE(*,*) "Current path OK? (Y or N): ",path
	DO
		READ(*,*) choice
		IF ((choice == 'y').OR.(choice == 'Y')) THEN
			EXIT
		ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
			WRITE(*,*) "Enter new path: "
			READ(*,'(a)') path
			path = TRIM(path)
			EXIT
		ELSE
			WRITE(*,*) 'Y or N please.'
		ENDIF
	END DO
END IF

DO
	WRITE(*,*) 'Enter file name (without extention -- .dat added): '
	READ(*,*) name
	filename = TRIM(path)//TRIM(name)//'.dat'
	OPEN(UNIT = 10, FILE = filename, STATUS = 'NEW',IOSTAT = iostatus)
	IF (iostatus > 0) THEN
		WRITE(*,*) "File already exists: Overwrite (o), append (a), new file name (n)?"
		READ(*,*) choice
		DO
			IF ((choice == 'o').OR.(choice == 'O')) THEN
				CLOSE(UNIT = 10)
				OPEN(UNIT = 10, FILE = filename, STATUS = 'REPLACE')
				EXIT
			ELSEIF ((choice == 'a').OR.(choice == 'A')) THEN
				CLOSE(UNIT = 10)
				OPEN(UNIT = 10,FILE = filename, STATUS = 'OLD', POSITION = 'APPEND')
				EXIT
			ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
				CLOSE(UNIT = 10)
				EXIT
			ELSE
				WRITE(*,*) 'Incorrect choice.  Please select again.'
			END IF
		END DO
		IF ((choice.NE.'n').AND.(choice.NE.'N')) EXIT
	ELSE
		EXIT
	END IF
END DO

DO i = 1, nf
	WRITE(10,1000) i,pevals(1,i),pevals(2,i),jvals(i)
END DO

CLOSE(UNIT = 10)

WRITE(*,*) ''
WRITE(*,*) 'Data written to:',filename

END SUBROUTINE J_WritePair

!====================================================================

SUBROUTINE PairLog(pairTEST)

USE phf_vals
USE system_parameters
USE spstate
USE psis

IMPLICIT NONE

LOGICAL, INTENT(OUT) :: pairTEST
INTEGER :: pairOP
INTEGER :: i,oddP,evenP

evenP = 0
oddP = 0
DO i = 1, nsps
	pairOP = (-1)**(spsqn(i)%l)
	IF (pairOP == 1) THEN
		evenP = evenP + 1
	ELSE
		oddP = oddP + 1
	END IF
END DO
WRITE(*,*)''
WRITE(*,*)'Checking Parity:'
IF (oddP == 0) THEN
	WRITE(*,*) ''
	WRITE(*,*) 'All Positive - ignoring'
	WRITE(*,*) ''
	pairTEST = .FALSE.
	npair = 1
ELSEIF (evenP == 0) THEN
	WRITE(*,*) ''
	WRITE(*,*) 'All Negative - ignoring'
	WRITE(*,*) ''
	pairTEST = .FALSE.
	npair = 1
ELSE
	WRITE(*,*) ''
	WRITE(*,*) 'Mixed: Including +/- parity'
	WRITE(*,*) ''
	pairTEST = .TRUE.
	npair = 2
END IF

END SUBROUTINE PairLog
