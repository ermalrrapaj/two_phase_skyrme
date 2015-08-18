	SUBROUTINE fdjac(x,fvec,df,usrfunc,bad)
	USE nr, ONLY: vecfunc
  USE nrtype; USE nrutil, ONLY : assert_eq
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: fvec
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: df
	LOGICAL, INTENT(OUT) :: bad
  PROCEDURE(vecfunc) :: usrfunc
	REAL(SP), PARAMETER :: EPS=1.0d-4
	INTEGER(I4B) :: j,n
	REAL(SP), DIMENSION(size(x)) :: xsav,xph,h,xpu,xpl
  
  bad = .false. 
	n=assert_eq(size(x),size(fvec),size(df,1),size(df,2),'fdjac')
	xsav=x
	h=EPS*abs(xsav)
	where (h == 0.0) h=EPS
	xph=xsav+h
	h=xph-xsav
  do j=1,n
    xpu = x 
    xpl = x
    if (abs(x(j))>EPS) then 
      xpu(j) = xpu(j)*(1.d0 + EPS) 
      xpl(j) = xpl(j)*(1.d0 - EPS) 
    else 
		  xpu(j) = xpu(j) + EPS 
      xpl(j) = xpl(j) - EPS
    endif    
    df(:,j)=(usrfunc(xpu)-usrfunc(xpl))/(xpu(j)-xpl(j)) 
    !x(j)=xph(j)
		!df(:,j)=(usrfunc(x)-fvec(:))/h(j)
		!x(j)=xsav(j)
	end do
  
  if (ANY(df.ne.df)) then
    bad = .true.
  !  print *,'bad Jacobian' 
  !  print *,df
  !  stop
  endif 
  
  END SUBROUTINE fdjac
