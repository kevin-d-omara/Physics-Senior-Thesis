!====================================================================
!
!  revised CWJ June 2015;
!  more logical output; sorts on energies
!  also prints out probabilities
! 
!  INPUT:
!     np      :  master declared dimension of arrays
!     nstates : = # of states
!     jval(:) : array of jvalues
!     eval(:) : array of energies
!     pval(:) : array of parities (+,-1)
!     parityflag : logical flag if true then possible to have different parities
!

SUBROUTINE J_WriteResults(tol,np,nf,jvals,pevals,nlist,jlist,problist,hamlist,parityflag)

use errortests
implicit none

real,    intent(in) :: tol
integer (kind=8), intent(in) :: np, nf
real    (kind=8), intent(in) :: jvals(np), pevals(2,np)
integer, intent(in) :: nlist
real    (kind=8), intent(in) :: jlist(nlist),problist(2,nlist),hamlist(2,nlist)
logical, intent(in) :: parityflag

real(kind=8) :: probsum,hamSum
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

!IF (path.EQ."undefined") THEN
!	WRITE(*,*) "Enter file path (i.e. ~/Desktop/path/)"
!	READ (*,'(a)') path
	!path = "'"//TRIM(path)//"'"
	!ELSE
	!WRITE(*,*) "Current path OK? (Y or N): ",path
	!DO
	!	READ(*,*) choice
	!	IF ((choice == 'y').OR.(choice == 'Y')) THEN
	!		EXIT
	!	ELSEIF ((choice == 'n').OR.(choice == 'N')) THEN
	!		WRITE(*,*) "Enter new path: "
	!		READ(*,'(a)') path
	!		path = TRIM(path)
	!		EXIT
	!	ELSE
	!		WRITE(*,*) 'Y or N please.'
	!	ENDIF
	!END DO
	!END IF

DO
	WRITE(*,*) 'Enter file name (without extention -- .dat added): '
	READ(*,*) name
	filename = TRIM(name)//'.dat'
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
write(10,*)' Easy Projection Results '
write(10,*)' (cutoff criterion = ',tol,')'
write(10,*)' '
if(parityflag)then
	write(10,*)' State    E(+ parity)  E(-parity)   J'
	write(10,*)' ------------------------------------'
	DO i = 1, nf
		WRITE(10,1001) i,pevals(1,i),pevals(2,i),jvals(i)		
	END DO
	1001 format(i3,2x,2F13.5,2x,f4.1)		

else
	write(10,*)' State    E       J'
	write(10,*)' ------------------------------------'
	DO i = 1, nf
		WRITE(10,1002) i,pevals(1,i),jvals(i)		
	END DO
1002 format(i3,2x,F12.5,2x,f4.1)	
end if	
!DO i = 1, nf
!	WRITE(10,1000) i,pevals(1,i),pevals(2,i),jvals(i)
!END DO
probsum = 0.d0 
write(10,*)' '
write(10,*)' Fraction in original HF state: Norm'
write(10,*)' '
if(parityflag)then
	write(10,*)' J    frac(+)   frac(-)'
	write(10,*)'-----------------------'
	do i = 1,nlist
		write(10,2001)jlist(i),problist(1,i),problist(2,i)
		probsum = probsum+problist(1,i)+problist(2,i)
	end do
	2001 format(f4.1,2f10.6)
else
	write(10,*)' J    frac'
	write(10,*)'-----------------------'
	do i = 1,nlist
		write(10,2001)jlist(i),problist(1,i)
		probsum = probsum+problist(1,i)
	end do
end if
write(10,*)' Total of HF state = ',probsum

hamsum = 0.d0
write(10,*)' '
write(10,*)' Fraction in original HF state: Hamiltonian'
write(10,*)' '
if(parityflag)then
	write(10,*)' J    frac(+)   frac(-)'
	write(10,*)'-----------------------'
	do i = 1,nlist
		write(10,2002)jlist(i),hamlist(1,i),hamlist(2,i)
		hamsum = hamsum+hamlist(1,i)+hamlist(2,i)
	end do
	2002 format(f4.1,2(1X,f10.6))
else
	write(10,*)' J    frac'
	write(10,*)'-----------------------'
	do i = 1,nlist
		write(10,2002)jlist(i),hamlist(1,i)
		hamsum = hamsum+hamlist(1,i)
	end do
end if
write(10,*)' Total of sum of trace(H) = ',hamSum

CLOSE(UNIT = 10)

WRITE(*,*) ''
WRITE(*,*) 'Data written to:',filename
return

END SUBROUTINE J_WriteResults
