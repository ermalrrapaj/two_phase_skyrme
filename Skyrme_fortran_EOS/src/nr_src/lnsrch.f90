	SUBROUTINE lnsrch(xold,fold,g,p,x,f,stpmax,check,func,usrfunc)
  USE nr_interface, ONLY: vecfunc,vecfunc_to_scalarfunc
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,vabs
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: xold,g
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
	REAL(SP), INTENT(IN) :: fold,stpmax
	REAL(SP), DIMENSION(:), INTENT(OUT) :: x
	REAL(SP), INTENT(OUT) :: f
	LOGICAL(LGT), INTENT(OUT) :: check
	!INTERFACE
	!	FUNCTION func(x,usrfunc)
	!	USE nrtype
  !  USE nr
	!	IMPLICIT NONE
	!	REAL(SP) :: func
	!	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	!  PROCEDURE(vecfunc) :: usrfunc
	!	END FUNCTION func
	!END INTERFACE
	PROCEDURE(vecfunc_to_scalarfunc) :: func
  PROCEDURE(vecfunc) :: usrfunc
  !REAL(SP), PARAMETER :: ALF=1.0e-4_sp,TOLX=epsilon(x)
	REAL(SP), PARAMETER :: ALF=1.0d-4,TOLX=epsilon(x)
	INTEGER(I4B) :: ndum,i
	REAL(SP) :: a,alam,alam2,alamin,b,disc,f2,pabs,rhs1,rhs2,slope,tmplam
	
  ndum=assert_eq(size(g),size(p),size(x),size(xold),'lnsrch')
	check=.false.
	pabs=vabs(p(:))
	if (pabs > stpmax) p(:)=p(:)*stpmax/pabs
	slope=dot_product(g,p)
	if (slope >= 0.0) call nrerror('roundoff problem in lnsrch')
	alamin=TOLX/maxval(abs(p(:))/max(abs(xold(:)),1.0_sp))
	alam=1.0
	! Check that we are not changing x too much over the step 
  !do i=1,SIZE(x) 
  !  if (abs(p(i))*alam/2.d0>abs(x(i))) then
  !    alam = abs(x(i))*2.d0/abs(p(i)) 
  !  endif 
  !enddo  
   
  do
		x(:)=xold(:)+alam*p(:)
		f=func(x,usrfunc)
		if (alam < alamin) then
			x(:)=xold(:)
			check=.true.
			RETURN
		else if (f <= fold+ALF*alam*slope) then
			RETURN
		else
			if (alam == 1.0) then
				tmplam=-slope/(2.0_sp*(f-fold-slope))
			else
				rhs1=f-fold-alam*slope
				rhs2=f2-fold-alam2*slope
				a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
				b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/&
					(alam-alam2)
				if (a == 0.0) then
					tmplam=-slope/(2.0_sp*b)
				else
					disc=b*b-3.0_sp*a*slope
					if (disc < 0.0) then
						tmplam=0.5_sp*alam
					else if (b <= 0.0) then
						tmplam=(-b+sqrt(disc))/(3.0_sp*a)
					else
						tmplam=-slope/(b+sqrt(disc))
					end if
				end if
				if (tmplam > 0.5_sp*alam) tmplam=0.5_sp*alam
			end if
		end if
		alam2=alam
		f2=f
		alam=max(tmplam,0.1_sp*alam)
	end do
	END SUBROUTINE lnsrch
