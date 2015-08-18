module eos_com_mod
  implicit none
  
  type eos_com
    real(8) :: nb  = -1.d0 
    real(8) :: np  = -1.d0 
    real(8) :: nn  = -1.d0 
    real(8) :: xp  = -1.d0 
    real(8) :: T   = -1.d0
    real(8) :: ee  = 0.d0
    real(8) :: pp  = 0.d0
    real(8) :: ss  = 0.d0 
    real(8) :: mun = 0.d0   
    real(8) :: mup = 0.d0   
    real(8) :: mue = 0.d0  
    real(8), dimension(30) :: xx = 0.d0 
    logical :: densities_set = .false. 
    logical :: cpotentials_set = .false. 
    logical :: full_set = .false.
    ! Mixed phase description 
    logical :: mixed_phase = .false.
    real(8) :: u = 0.d0  
    real(8) :: npl = 0.d0, npu = 0.d0
    real(8) :: nnl = 0.d0, nnu = 0.d0
  end type eos_com

contains

  function eos_input(nb,xp,T) result(eos) 
    real(8), intent(in) :: nb,xp
    real(8), intent(in),optional :: T 
    type(eos_com) :: eos
    eos%densities_set = .true. 
    eos%nb = nb 
    eos%xp = xp
    eos%np = nb*xp 
    eos%nn = nb*(1.d0-xp) 
    eos%T  = -1.d0
    if (PRESENT(T)) eos%T  = T  
  end function eos_input
  
  function eos_input_from_mixed(np,mun,T) result(eos) 
    real(8), intent(in) :: np,mun
    real(8), intent(in), optional :: T  
    type(eos_com) :: eos 
    eos%mun = mun 
    eos%np = np
    eos%T = -1.d0
    if (PRESENT(T)) eos%T = T
  end function eos_input_from_mixed 
   
  function eos_input_from_mu(mup,mun,T) result(eos) 
    real(8), intent(in) :: mup,mun
    real(8), intent(in), optional :: T  
    type(eos_com) :: eos 
    eos%cpotentials_set = .true. 
    eos%mun = mun 
    eos%mup = mup 
    eos%T = -1.d0
    if (PRESENT(T)) eos%T = T
  end function eos_input_from_mu 

  function eos_input_from_densities(np,nn,T) result(eos) 
    real(8), intent(in) :: np,nn
    real(8), intent(in), optional :: T  
    type(eos_com) :: eos 
    eos%densities_set = .true.
    eos%nn = nn 
    eos%np = np
    eos%nb = np + nn 
    eos%xp = np/(np + nn + 1.d-50)  
    eos%T = -1.d0
    if (PRESENT(T)) eos%T = T
  end function eos_input_from_densities 

end module eos_com_mod

