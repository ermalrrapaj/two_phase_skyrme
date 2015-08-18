module nr_interface
  
  INTERFACE
  	FUNCTION vecfunc(x,output_call) RESULT(Y)
  	USE nrtype
  	IMPLICIT NONE
  	REAL(SP), DIMENSION(:), INTENT(IN) :: x
  	LOGICAL, OPTIONAL :: output_call
    REAL(SP), DIMENSION(size(x)) :: Y
  	END FUNCTION vecfunc
  END INTERFACE
  
  INTERFACE
  	FUNCTION vecfunc_to_scalarfunc(x,usrfunc) RESULT(Y)
  	USE nrtype
  	IMPLICIT NONE
  	REAL(SP), DIMENSION(:), INTENT(IN) :: x
    PROCEDURE(vecfunc) :: usrfunc 
  	REAL(SP) :: Y
  	END FUNCTION vecfunc_to_scalarfunc
  END INTERFACE

end module nr_interface

