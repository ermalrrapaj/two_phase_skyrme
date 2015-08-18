program LStest 
  use eos_com_mod
  use LS_construction 
  
  type(eos_com) :: eosin 
  type(eos_com) :: eosout 
  
  eosin = eos_input(0.01d0,0.4d0,1.d0) 
  eosout = get_LS_state(eosin,.true.)
   
end program
