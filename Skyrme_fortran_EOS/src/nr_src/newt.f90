	SUBROUTINE newt(x,check,usrfunc)
	USE nrtype; USE nrutil, ONLY : nrerror,vabs
	USE nr, ONLY : vecfunc,fdjac,lnsrch,lubksb,ludcmp
	USE fminln
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: x
	LOGICAL(LGT), INTENT(OUT) :: check
	INTEGER(I4B), PARAMETER :: MAXITS=200
	REAL(SP), PARAMETER :: TOLF=1.0e-10_sp,TOLMIN=1.0e-10_sp,TOLX=epsilon(x),&
		STPMX=100.0
	INTEGER(I4B) :: its
	INTEGER(I4B), DIMENSION(size(x)) :: indx
	REAL(SP) :: d,f,fold,stpmax
	REAL(SP), DIMENSION(size(x)) :: g,p,xold
	REAL(SP), DIMENSION(size(x)), TARGET :: fvec
	REAL(SP), DIMENSION(size(x),size(x)) :: fjac
	PROCEDURE(vecfunc) :: usrfunc
  LOGICAL :: bad
	fmin_fvecp=>fvec
	f=fmin(x,usrfunc)
	if (maxval(abs(fvec(:))) < 0.01_sp*TOLF) then
		check=.false.
		RETURN
	end if
	stpmax=STPMX*max(vabs(x(:)),real(size(x),sp))
	do its=1,MAXITS
		call fdjac(x,fvec,fjac,usrfunc,bad)
		if (bad) then
      check = .true. 
      return 
    endif   
    g(:)=matmul(fvec(:),fjac(:,:))
		xold(:)=x(:)
		fold=f
		p(:)=-fvec(:)
		call ludcmp(fjac,indx,d)
		call lubksb(fjac,indx,p)
		call lnsrch(xold,fold,g,p,x,f,stpmax,check,fmin,usrfunc)
		if (maxval(abs(fvec(:))) < TOLF) then
      check=.false.
			RETURN
		end if
		if (check .and. abs(f)<TOLMIN) then
      check=(maxval(abs(g(:))*max(abs(x(:)),1.0_sp) / &
				max(f,0.5_sp*size(x))) < TOLMIN)
			RETURN
		end if
    if (maxval(abs(x(:)-xold(:))/max(abs(x(:)),1.0_sp)) < TOLX) &
			RETURN
	end do
	!print *,g
  check = .true.
  call nrerror('MAXITS exceeded in newt')
	END SUBROUTINE newt
