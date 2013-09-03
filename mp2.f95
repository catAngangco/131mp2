PROGRAM MP2
IMPLICIT NONE

character*40 FILENAME
character*2 answer, yes,no
integer :: lu, count, n, n1, itmax, irow, jcol, stat
double precision,dimension(:,:),allocatable :: A, I, O
double precision :: det,eg, en, cond1,cond2,cond8,condf

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
	allocate(O(n,n))

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

	print *,"A:", A

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
	O = A
	!inverse
	det = 1.0
	A = gaussjordan(A,det)
	WRITE(lu,*) "Inverse:"
	WRITE(lu,*) A
	print *,"I", A
	WRITE(lu,*) "Determinant:",det
	print *,"Det:",det

	!cond2
	cond2 = norm2(O)
	print *,"cond2:",cond2
	WRITE(lu,*) "Cond2:",cond2

	O=abs(O)
	A=abs(A)
	!cond1
	cond1 = norm1(O) * norm1(A)
	print *,"cond1:",cond1
	WRITE(lu,*) "Cond1:",cond1

	!cond8
	cond8 = norm8(O) * norm8(A)
	print *,"cond8:",cond8
	WRITE(lu,*) "Cond8:",cond8

	!condf
	condf = normf(O**2) * normf(A**2)
	print *,"condf:",condf
	WRITE(lu,*) "CondF:",condf

	deallocate(A)
	deallocate(O)
	1335 return
end subroutine matrix

double precision function norm2(A)
	double precision, dimension(n,n) :: A
	double precision, dimension(2) :: eval
	eval = getEigenvalues(A)
	print *,"E",eval
	norm2 = sqrt(eval(1)/eval(2))
end function
double precision function norm1(A)
	double precision, dimension(n,n) :: A
	double precision, dimension(n) :: colsums
	colsums = sum(A,DIM = 1)
	norm1 = maxval(colsums)
end function

double precision function norm8(A)
	double precision, dimension(n,n) :: A
	double precision, dimension(n) :: rowsums
	rowsums = sum(A,DIM = 2)
	norm8 = maxval(rowsums)
end function

double precision function normf(A)
	double precision, dimension(n,n) :: A
	normf = sqrt(sum(A))
end function

integer function pivotOf(A,k)
	double precision, dimension(n) :: r
	double precision, dimension(n,n) :: A
	integer :: i,k

	do i = 1,n
		r(i) = abs(A(k,i))
	end do
	pivotOf = maxloc(r,1)
end function

function gaussjordan(A,det)
	double precision, dimension(n,n) :: gaussjordan, A
	integer, dimension(n) :: p
	double precision, dimension(n) :: t
	integer :: i,j,k,loc
	double precision :: aik,temp,pivot, amax,det
	do i = 1,n
		p(i) = i
	end do
	do k = 1,n
		loc = pivotOf(A,k)
		amax = A(k,loc)
		if (abs(amax).lt.eg) then
			goto 10
		else
			det = det * amax
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
end function

double precision function tr(T,trLen)
integer :: i,trLen
double precision, dimension(trLen,trlen) :: T
double precision :: q
q=0
do i=1,trLen
	q = q + T(i,i)
end do
tr = q
end function

function multiI(qu,nu)
integer :: nu,iter, iter2
integer, dimension(nu,nu) :: I, multiI
double precision :: qu
do iter =1, nu
	do iter2=1,nu
		if (iter.eq.iter2) then
			I(iter,iter2) = qu
		else
			I(iter,iter2) = 0
		endif
	end do
end do
multiI = I
end function

function matrixSub(Ai,Bi,ni)
integer :: i,j,ni
double precision, dimension(ni,ni) :: matrixSub,Ai,Bi,Cii
do i = 1, ni
	do j =1,ni
		Cii(i,j) = Ai(i,j) - Bi(i,j)
	end do
end do
matrixSub = Cii
end function

function getBk(Ag,Bg,qg,ng)
integer :: ng
double precision :: qg
double precision, dimension(ng,ng):: Ag,Bg,Cig,Dg,Iqg, getBk
Iqg = multiI(qg,ng)
Cig = matrixSub(Bg,Iqg,ng)
Dg = matmul(Ag,Cig)
getBk = Dg
end function

double precision function f(Y,coefficients)
double precision :: Y
integer :: i
integer, dimension(n1) :: rev,coefficients
f=0.0
do i=0,n
	rev(n1-i) = i
end do
do i=1,n1
	!print *,i,f,"+",coefficients(i),"x" ,Y,"^",rev(i)
	f = f + (1.0*coefficients(i)*(Y**rev(i)))
	!print *,"=",f
end do

end function

function graeffe(Af)
integer :: Af(:)
double precision, dimension(size(AF)) :: T, Bf
double precision, dimension(size(AF)-1) :: groots, X, graeffe
integer:: kf,mf, switch, i, jf, lf, sign
double precision :: low, curr, fx, fy
double precision :: high
print *,""
print *,"GRSM:"
T = 0.0
Bf = 0.0
X = 0.0
groots = 0.0
low = 10D-150
high = 10D+150
kf = 0
mf = 1
Bf = Af
print *,Bf
switch = 1
10      continue
		kf = kf + 1
		do i = 2, size(AF)
			jf = i-1
			lf = i +1
			sign = -1
			T(i) = Bf(i)**2
			if ((jf.ge.1).and.(lf.le.size(AF))) then
				goto 30
			else
				goto 40
			endif
			30      continue
					!print *,i,T(i),Bf(jf),Bf(lf)
					T(i) = T(i) + (sign*2*Bf(jf)*Bf(lf))
					!print *,T(i)
					jf = jf-1
					lf = lf+1
					sign = -sign
					if ((jf.ge.1).and.(lf.le.size(AF))) then
						goto 30
					endif
					curr = ABS(T(i))
					if (curr.gt.high) then
						switch = 0
					endif
					if (curr.lt.low) then
						switch = 0
					endif
			40		continue
		end do
		mf = 2**kf
		do i=2,n1
			Bf(i) = T(i)
		end do
		print *,Bf
		if (switch.eq.1) then
			goto 10
		else
			goto 20
		endif
20      continue
!Calulate roots
do i=2,size(AF)
	X(i-1) = (ABS(Bf(i)/Bf(i-1)))**(1.0/mf)
end do
!print *,"X"
!print *,X
! Determine sign
do i=1,n
	80      continue
			!print *,"loop"
			fx = f(X(i),Af)
			if (abs(fx).lt.en) then
				!print *,"root:",X(i)
				groots(i) = X(i)
			elseif (X(i).gt.0) then
				X(i) = -X(i)
				goto 80
			else
				!print *,"not root:",X(i)
			endif
end do
if (groots(1).eq.0) then
	do i=1,n
		fx = f(X(i),Af)
		fy = f(-X(i),Af)
		if (abs(fx).lt.abs(fy)) then
			groots(i) = X(i)
		else
			groots(i) = -X(i)
		endif
	end do
endif
graeffe = groots
end function

double precision function df(dir, coefficients)
double precision :: dir
integer, dimension(n1) :: rev, coefficients
integer :: i
df = 0.0
do i=0,n
	rev(n1-i) = i
end do
do i=1,n1
	df = df + (coefficients(i)*rev(i)*(dir**(rev(i)-1)))
end do
end function

double precision function newton(w,coefficients)
integer :: currItr
double precision :: dfw,fw, w
integer, dimension(n1) :: coefficients
currItr = 0
fw = f(w,coefficients)
print *,""
print *,currItr,w,fw,"Initial estimate"
do currItr=1,itmax
  	!print *,"iter",w
	dfw = df(w,coefficients)
    !print *,dfw
	if (dfw.eq.0) then 
		newton = w
		goto 50
	endif
	w = w - (fw/dfw)
	fw = f(w,coefficients)
	print *,currItr,w,fw
	if (abs(fw) < en) then
		newton = w
		!print *,"e"
		goto 50
	endif
end do
newton = w
50      continue
end function

function getEigenvalues(A)
	double precision, dimension(n,n) :: A, B
	double precision, dimension(2) :: getEigenvalues, E
	double precision, dimension(n) :: roots, nroots
	integer, dimension(n1) :: coefficients
	integer :: k
	double precision :: q

	! c1
	coefficients(1) = 1
	q = tr(A,n)
	coefficients(2) = -q
	! c2
	B=getBk(A,A,q,n)
	q=tr(B,n)
	q=q/2
	coefficients(3) = -q
    !c3
    do k=3,n
        B=getBk(A,B,q,n)
		q=tr(B,n)
		q=q/k
		coefficients(k+1) = -q
	end do
	print *,""
	print *,"COEFFICIENTS"
	print *,coefficients
	!call printPoly(coefficients)
	!GRSM
	roots = graeffe(coefficients)
	print *,"Graeffe roots"
	print *,roots
	print *,""
	print *,"NEWTON"
	do irow=1,n
		nroots(irow) = newton(roots(irow),coefficients)
	end do
	print *,"Newton roots"
	print *,nroots

	getEigenvalues(1) = maxval(nroots)
	getEigenvalues(2) = minval(nroots)
end function
END PROGRAM





