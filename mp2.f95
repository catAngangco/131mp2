PROGRAM MP2
IMPLICIT NONE

character*40 FILENAME
character*2 answer, yes,no
integer :: lu, count, n, n1, itmax, irow, jcol, stat
double precision,dimension(:,:),allocatable :: A, I
double precision :: eg, en

WRITE(*,*) "Enter input file pathname:"
READ(*,*) FILENAME
lu=9
open(lu,file=FILENAME)
no = "N"
do
	1336 call matrix
	1337 WRITE(*,*) "Read the next entry (Y/N)"
	READ(*,*) answer
	if (answer.eq.'n'.or.answer.eq.'N') exit
	if (answer.ne.'y'.and.answer.ne.'Y') goto 1337
	goto 1336
end do
70 continue
close(9)
close(30)

contains
subroutine matrix
	READ(9,*, IOSTAT = stat) n
	IF (stat < 0) THEN
    print*, "END OF FILE"
    GOTO 1335
    END IF
	n1 = n + 1

	allocate(A(n,n))
	allocate(I(n,n))

	READ(9,*, IOSTAT = stat) A
	IF (stat < 0) THEN
    print*, "END OF FILE"
    GOTO 1335
    END IF
    
	READ(9,*, IOSTAT = stat) eg
	IF (stat < 0)  THEN
    print*, "END OF FILE"
    GOTO 1335
    END IF

    READ(9,*, IOSTAT = stat) en
	IF (stat < 0)  THEN
    print*, "END OF FILE"
    GOTO 1335
    END IF
    
	READ(9,*, IOSTAT = stat) itmax
	IF (stat < 0)  THEN
    print*, "END OF FILE"
    GOTO 1335
    END IF

	print *,"A"
	do irow =1,n
		print*, (A(irow,jcol), jcol=1,n)
	end do

	lu=30
	open(lu,file = "outputData.txt")
	WRITE(lu,*) "Entry:"
	WRITE(lu,*) "A:"
	do irow =1,n
		WRITE(lu,*) (A(irow,jcol), jcol=1,n)
	end do
	WRITE(lu,*) "Eg:", eg
	WRITE(lu,*) "En:", en
	WRITE(lu,*) "itmax:",itmax

	I = gaussjordan(A)

	print *,"I:",I

	deallocate(A)
	deallocate(I)
	1335 return
end subroutine matrix

integer function pivotOf(A,k)
	double precision, dimension(n) :: r
	double precision, dimension(n,n) :: A
	integer :: i,k

	do i = 1,n
		r(i) = abs(A(k,i))
	end do
	pivotOf = maxloc(r,1)
end

function gaussjordan(A)
	double precision, dimension(n,n) :: gaussjordan, A
	integer, dimension(n) :: p
	double precision, dimension(n) :: t
	integer :: i,j,k,loc
	double precision :: aik,temp,pivot, amax
	do i = 1,n
		p(i) = i
	end do

	do k = 1,n
		loc = pivotOf(A,k)
		amax = A(k,loc)
		if (abs(amax).lt.eg) then
			goto 10
		endif
		!swap rows
		if (loc.ne.k) then
			do i = 1,n
				temp=A(i,k)
				A(i,k)=A(i,loc)
				A(i,loc)=temp
			end do
			temp = p(k)
			p(k) = p(loc)
			p(loc) = temp
		endif
		!normalize
		do j=1,n
			A(j,k)=A(j,k)/amax
		end do
		A(k,k)=1.0/amax
		!reduce
		do i=1,n
			if(i.ne.k) then
				aik = A(k,i)
				do j=1,n
					A(j,i)=A(j,i)-(aik*A(j,k))
				end do
				A(k,i) = -aik/amax
			endif
		end do
	end do
	!unscramble

	do i=1,n
		do j=1,n
			t(p(j)) = A(j,i)
		end do
		do j=1,n
			A(j,i) = t(j)
		end do
	end do
	10      continue
	gaussjordan=A
end

END PROGRAM





