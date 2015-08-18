program tabular_test
  use tabular_EOS 
  use eos_com_mod 
  implicit none 

  type(eos_com) :: eos 
  real(8) :: np,nn,mun,T 
  integer :: i  
  call read_table("./eostable_Skyrme.h5") 
  T = T_arr(10) 
  np = 1.d-3
  do i=1,501
    nn = 10.d0**((LOG10(1.d0)-LOG10(1.d-10))*DBLE(i-1)/500.d0 + LOG10(1.d-10)) 
    eos = get_bulk_state(eos_input_from_densities(np,nn,T)) 
    print *,eos%nn,eos%mup,eos%mun,eos%pp 
  enddo   
  
  eos = get_bulk_state_from_mixed(eos_input_from_mixed(1.d-4,-2.d0*T,1.d0/197.3d0))  
   

end program 
