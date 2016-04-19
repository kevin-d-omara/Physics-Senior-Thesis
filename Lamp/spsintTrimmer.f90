! ##################################################################################################
! PURPOSE:
!	Automate task of trimming .sps and .int files.  Simply input the filename and the highest
!	Harmonic Oscillator Shell found in the .sps file.  The program will create all trimmed versions
!	of this file, down to the highest H.O. Shell of 0.  May also be used to trim .int files
!	simultaneously.
! ##################################################################################################

MODULE globalVar

	IMPLICIT NONE
	INTEGER :: unitNumOld, unitNumNew, maxHO
	INTEGER :: unitNumOldINT, unitNumNewINT, numHeaderOld, numHeaderNew
	INTEGER :: numOrbitsNew, numOrbitsOld
	INTEGER :: numLinesNew, numLinesOld
	CHARACTER(LEN=20) :: filenameSPS, filenameOld, filenameNew
	CHARACTER(LEN=20) :: filenameINT, filenameOldINT, filenameNewINT
	CHARACTER(LEN=20) :: tempString
	CHARACTER(LEN=60) :: tempStringLong
	CHARACTER(LEN=150) :: tempStringLongLong
	CHARACTER(LEN=1) :: choice
	LOGICAL :: trimINT

END MODULE globalVar

! -------------------------------------------------------------

PROGRAM spsintTrimmer
	IMPLICIT NONE

 CALL mainMenu()
 CALL mainBody()
 WRITE(*,*) 'Completed'

END PROGRAM spsintTrimmer



! ##################################################################################################
!
! Use enters name of file to be trimmed and the highest harmonic oscillator (HO) shell contained
!   within.  Optionally enters information about the .int file.
!
! ##################################################################################################
SUBROUTINE mainMenu()
	USE globalVar

WRITE(*,*) ' Welcome to the .sps & .int file trimmer.'
WRITE(*,*) "           by Kevin O'Mara"
WRITE(*,*) '      Version 0.4 : July 1, 2015'
WRITE(*,*) ''
WRITE(*,*) ' * * * * * * * * * * * * * * * * * * * *'
WRITE(*,*) ''

WRITE(*,*) "Enter the name of the file to be trimmed (.sps)"						! SPS filename
READ(*,*) filenameSPS
OPEN(unit=1, file=(TRIM(ADJUSTL(filenameSPS)) // '.sps'), status='OLD')
!IF ERROR, LOOP

WRITE(*,*) "What do you want the new files to be named? (.sps)"						! SPS filename New
WRITE(*,*) "These will be appended with '#', the number of max orbitals."
READ(*,*) filenameSPS

filenameOld = filenameSPS
unitNumOld = 1

WRITE(*,*) "Enter the highest harmonic oscillator shell this file contains (N)."	! enter maxHO
READ(*,*) maxHO

WRITE(*,*) 'Do you also want to trim the accompanying .int file? (Y/N)'				! trim .int ?
DO
    READ(*,*) choice
    IF ((choice == 'y').OR.(choice == 'Y')) THEN
	    trimINT = .TRUE.
        EXIT
    ELSE IF ((choice == 'n').OR.(choice == 'N')) THEN
	    trimINT = .FALSE.
        EXIT
    ELSE
	    WRITE(*,*) 'Y or N please.'
    END IF
END DO

IF (trimINT.EQV..TRUE.) THEN
	WRITE(*,*) "Enter the name of the file to be trimmed (.int)"					!INT filename
	READ(*,*) filenameINT
	OPEN(unit=2, file=(TRIM(ADJUSTL(filenameINT)) // '.int'), status='OLD')
	!IF ERROR, LOOP

	WRITE(*,*) "What do you want the new files to be named? (.int)"					!INT filename New
	WRITE(*,*) "These will be appended with '#', the number of max orbitals."
	READ(*,*) filenameINT

	filenameOldINT = filenameINT
	unitNumOldINT = 2

	WRITE(*,*) "Enter the number of header lines in the .int file."					! enter numHeaderOld
	READ(*,*) numHeaderOld
END IF


RETURN
END SUBROUTINE mainMenu



! ##################################################################################################
!
!	Builds new .sps files and optionally .int files.
!
! ##################################################################################################
SUBROUTINE mainBody()
	USE globalVar
	IMPLICIT NONE
	INTEGER :: i

DO i = maxHO,1,-1			!Work backwards from current HO shell.
	maxHO = maxHO - 1			!Start with next highest harmonic oscillator (HO) shell.

	CALL createNewSPS()
	IF (trimINT.EQV..TRUE.) CALL createNewINT()

	CALL findNumOrbits()
	IF (trimINT.EQV..TRUE.) CALL findNumLines()

	CALL writeNewSPS()
	IF (trimINT.EQV..TRUE.) CALL writeNewINT()

	CALL updateFileIndex()
END DO

RETURN
END SUBROUTINE mainBody



! ##################################################################################################
!
! Builds new .sps file, labled "filename#.sps" (i.e. maxorb6.sps). 
!
!	INPUT   :: maxHO, filenameSPS
!
!	RETURNS :: unitNumNew, filenameNew
!
! ##################################################################################################
SUBROUTINE createNewSPS()
	USE globalVar

CALL integerTOstring(maxHo,tempString)
WRITE(*,*) 'Creating .sps file with maximum harmonic oscillator shell: ', TRIM(ADJUSTL(tempString))

unitNumNew = 1000 + maxHO										! Creates current working file.
filenameNew = TRIM(ADJUSTL(filenameSPS)) // TRIM(ADJUSTL(tempString))
OPEN(unit=unitNumNew, file=(TRIM(ADJUSTL(filenameNew)) // '.sps'), status='UNKNOWN')

RETURN
END SUBROUTINE createNewSPS



! ##################################################################################################
!
! Takes an integer from input -> returns string as output.
!
!	INPUT   :: 'num' (any integer value) & 'string' (any string of length 20)
!
!	RETURNS :: 'string' (usually 'tempString')
!
! ##################################################################################################
SUBROUTINE integerTostring(num,string)
	IMPLICIT NONE
	INTEGER, INTENT(in) :: num
	CHARACTER(LEN=20) :: string

WRITE(string,'(I20)') num

RETURN
END SUBROUTINE integerTostring



! ##################################################################################################
!
! Takes a real value from input -> returns string as output.
!
!	INPUT   :: 'num' (any real value) & 'string' (any string of length 20)
!
!	RETURNS :: 'string' (usually 'tempString')
!
! ##################################################################################################
SUBROUTINE realTostring(num,string)
	IMPLICIT NONE
	REAL, INTENT(in) :: num
	CHARACTER(LEN=20) :: string

WRITE(string,'(F12.6)') num

RETURN
END SUBROUTINE realTostring



! ##################################################################################################
!
! Builds new .int file, labled "filename#.int" (i.e. A12maxor8lam6.sps).
!									    ^						 ^
!	INPUT   :: maxHO, filenameINT
!
!	RETURNS :: unitNumNewINT, filenameNewINT
!
! ##################################################################################################
SUBROUTINE createNewINT()
	USE globalVar

CALL integerTOstring(maxHo,tempString)
WRITE(*,*) 'Creating .int file with maximum harmonic oscillator shell: ', TRIM(ADJUSTL(tempString))

unitNumNewINT = 2000 + maxHO										! Creates current working file.
filenameNewINT = TRIM(ADJUSTL(filenameINT)) // TRIM(ADJUSTL(tempString))
OPEN(unit=unitNumNewINT, file=(TRIM(ADJUSTL(filenameNewINT)) // '.int'), status='UNKNOWN')

RETURN
END SUBROUTINE createNewINT



! ##################################################################################################
!
! Finds the number of orbits for the new .sps file.
!
!	RETURNS :: numOrbitsNew
!
! ##################################################################################################
SUBROUTINE findNumOrbits()
	USE globalVar
	REAL :: a, b, c

numOrbitsNew = 0
READ(unitNumOld,*); READ(unitNumOld,*)			! Skip first two lines
DO												! Find new number of orbits
	READ(unitNumOld,*) a, b, c, numOrbitsOld
	IF (numOrbitsOld > maxHO) THEN
		EXIT
	ELSE
		numOrbitsNew = numOrbitsNew + 1
	END IF
END DO

REWIND(unitNumOld)
READ(unitNumOld,*); READ(unitNumOld,*)			! Skip first two lines
	
RETURN
END SUBROUTINE findNumOrbits



! ##################################################################################################
!
! Finds the number of lines (both total and header) for the new .int file.
!
!	RETURNS :: numLinesNew, numHeaderNew
!
! ##################################################################################################
SUBROUTINE findNumLines()
	USE globalVar
	INTEGER :: i, j
	INTEGER :: a, b, c, d

numHeaderNew = 0
j = 0
DO i = 1, numOrbitsNew						! Find new number of header lines
	j = j + 1
	IF (MOD(j,10)==0) THEN
		numHeaderNew=numHeaderNew+1
		j = 0
	END IF
END DO
IF (j>0) numHeaderNew = numHeaderNew+1		! Account for trailing header line

numLinesNew = 0
READ(unitNumOldINT,*) numLinesOld
DO i = 1, numHeaderOld-1					! Skip header lines
	READ(unitNumOldINT,*)
END DO

DO i = 1, (numLinesOld-numHeaderOld-1)		! Find new number of lines
	READ(unitNumOldINT,*) a, b, c, d
	IF ((a>numOrbitsNew).OR.(b>numOrbitsNew).OR.(c>numOrbitsNew).OR.(d>numOrbitsNew)) THEN
		CONTINUE
	ELSE
		numLinesNew = numLinesNew + 1
	END IF
END DO

REWIND(unitNumOldINT)
	
RETURN
END SUBROUTINE findNumLines



! ##################################################################################################
!
! Writes body of new .sps file.
!
! ##################################################################################################
SUBROUTINE writeNewSPS()
	USE globalVar
	IMPLICIT NONE
	INTEGER :: j

WRITE(unitNumNew,'(A)') 'iso'								! Writes first two lines: 'iso' and number of orbits.
CALL integerTostring(numOrbitsNew,tempString)
WRITE(unitNumNew,'(10X,A)') TRIM(ADJUSTL(tempString))

DO j = 1,numOrbitsNew										! Writes body of .sps file.
	READ(unitNumOld,'(A)') tempString
	WRITE(unitNumNew,'(A)') tempString
END DO

RETURN
END SUBROUTINE writeNewSPS



! ##################################################################################################
!
! Writes body of new .int file.
!
!	INPUT   :: numLinesOld
!
! ##################################################################################################
SUBROUTINE writeNewINT()
	USE globalVar
	IMPLICIT NONE
	INTEGER :: i, j
	INTEGER :: a, b, c, d
	REAL, DIMENSION(10) :: k

CALL writeHeader()

!READ(unitNumOldINT,*) numLinesOld			! Skip old header
DO i = 1, numHeaderOld-1
	READ(unitNumOldINT,*) 
END DO

DO j = 1, numLinesOld-numHeaderOld-1		! Only write lines with values <= numOrbitNew
	READ(unitNumOldINT,*) a, b, c, d
	IF ((a>numOrbitsNew).OR.(b>numOrbitsNew).OR.(c>numOrbitsNew).OR.(d>numOrbitsNew)) THEN
		CONTINUE
	ELSE
		BACKSPACE(unitNumOldINT)
		READ(unitNumOldINT,'(A)') tempStringLong
		WRITE(unitNumNewINT,'(A)') tempStringLong
	END IF
END DO

WRITE(unitNumNewINT,*) numLinesNew			! Write trailing line, contains numLinesNew

RETURN
END SUBROUTINE writeNewINT



! ##################################################################################################
!
! Writes header of new .int file.
!
!	RETURNS :: numLinesOld
!
! ##################################################################################################
SUBROUTINE writeHeader()
	USE globalVar
	IMPLICIT NONE
	INTEGER :: i, n, counter
	LOGICAL :: lineOne
	REAL, DIMENSION(10) :: f

CALL integerTOstring(numLinesNew,tempString)		! Prep for header writing
tempStringLongLong = TRIM(ADJUSTL(tempString))
lineOne = .TRUE.
counter = numOrbitsNew

DO
	IF (counter >= 10) THEN								! Determine number of values on next line
		n = 10
	ELSE IF (counter == 0 ) THEN
		EXIT
	ELSE
		n = counter
	END IF

	IF (lineOne.EQV..TRUE.) THEN						! Read values from unitNumNewINT into array f(10)
		READ(unitNumOldINT,*) numLinesOld, (f(i), i=1,n)
		lineOne = .FALSE.
	ELSE
		READ(unitNumOldINT,*) (f(i), i=1,n)
	END IF

	DO i = 1, n											! Append values onto tempStringLongLong
		CALL realTOstring(f(i),tempString)
		tempStringLongLong = TRIM(ADJUSTL(tempStringLongLong)) // '  ' // TRIM(ADJUSTL(tempString))
	END DO

	WRITE(unitNUmNewINT,'(2X,A)') TRIM(ADJUSTL(tempStringLongLong))		! Write to file
	tempStringLongLong = ''												! "Flush buffer"
	counter = counter - n
END DO


RETURN
END SUBROUTINE writeHeader



! ##################################################################################################
!
! Updates the unitNum of .sps & .int files to prepare for next round of trimmings.
!
!	INPUT   :: unitNumNew, unitNumNewINT, numHeaderNew
!
!	RETURNS :: unitNumOld, unitNumOldINT, numHeaderOld
!
! ##################################################################################################
SUBROUTINE updateFileIndex()
	USE globalVar

 CLOSE(unitNumOld); CLOSE(unitNumNew)
unitNumOld = unitNumNew								! Re-open the 'new' file (now it is old).
OPEN(unit=unitNumNew, file=(TRIM(ADJUSTL(filenameNew)) // '.sps'), status='OLD')

IF (trimINT.EQV..TRUE.) THEN
	CLOSE(unitNumOldINT); CLOSE(unitNumNewINT)
	unitNumOldINT = unitNumNewINT
	OPEN(unit=unitNumNewINT, file=(TRIM(ADJUSTL(filenameNewINT)) // '.int'), status='OLD')
	numHeaderOld = numHeaderNew
END IF

RETURN
END SUBROUTINE updateFileIndex
