module bisect_mod
  implicit none
 
contains 
  function bisect(xmin,xmax,check,usrfunc) result(xx)
    real(8) :: xx 
    real(8), intent(in)  :: xmin
    real(8), intent(in), optional :: xmax
    logical, intent(out) :: check
	  INTERFACE
	  	FUNCTION usrfunc(x) result(yy)
	  	IMPLICIT NONE
      real(8) :: yy
	  	real(8), INTENT(IN) :: x
	  	END FUNCTION usrfunc
	  END INTERFACE

    real(8) :: yl,ym,yu 
    real(8) :: xl,xm,xu
    integer :: i 

    check = .true.  
    xl = xmin 
    if (PRESENT(xmax)) then 
      xu = xmax 
    else 
      xu = xmin ! This will force NR iteration 
    endif  
    yl = usrfunc(xl)  
    yu = usrfunc(xu)
    
    if (yl*yu>0.d0) then
      !print *,'Bisect got a bad interval',yl,yu
      !return
      i = 1
      xm = xmin
    else  
      do i=1,50  
        xm = 0.5d0*(xu + xl) 
        ym = usrfunc(xm) 
        if (ym*yu>=0.d0) then
          yu = ym 
          xu = xm 
        else 
          yl = ym 
          xl = xm 
        endif 
        if (ABS(yu-yl)<1.d-3) exit     
      enddo  
    endif  
    
    ! Polish off with NR if bisection got close
    if (i<50) then 
      do i=1,20 
        yu = usrfunc(xm*(1.d0 + 1.d-6)) 
        yl = usrfunc(xm*(1.d0 - 1.d-6)) 
        ym = usrfunc(xm) 
        xm = xm - ym/(yu-yl)*2.d-6*xm
        if (ABS(ym)<1.d-12) exit   
      enddo  
      if (ABS(ym)>1.d-12) then ! NR failure condition
        xx = xm
        check = .true. 
        return 
      endif  
    else 
      !print *,'bisect probably failed',yl,yu,xl,xu 
      xx = xm
      return 
    endif 
    
    xx = xm 
    check = .false.  
    
  end function bisect
end module bisect_mod


