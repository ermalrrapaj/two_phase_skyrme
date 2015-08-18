MODULE fminln
	USE nr_interface, ONLY: vecfunc
  USE nrtype; USE nrutil, ONLY : nrerror
	REAL(SP), DIMENSION(:), POINTER :: fmin_fvecp
CONTAINS
!BL
	FUNCTION fmin(x,usrfunc)
	  IMPLICIT NONE
	  REAL(SP), DIMENSION(:), INTENT(IN) :: x
	  REAL(SP) :: fmin
    PROCEDURE(vecfunc) :: usrfunc
	  if (.not. associated(fmin_fvecp)) call &
	  	nrerror('fmin: problem with pointer for returned values')
	  fmin_fvecp=usrfunc(x,.true.)
	  fmin=0.5_sp*dot_product(fmin_fvecp,fmin_fvecp)
	END FUNCTION fmin
END MODULE fminln
