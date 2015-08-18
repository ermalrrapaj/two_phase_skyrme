program nse_test 
  use eos_nse_mod
  use eos_com_mod
   
  implicit none 
  
  type(eos_com) :: eos,eos_outside
  type(nucleus) :: nuclei(50) 
  integer :: i 
  real(8) :: nb,xp,T 
  
  nb = 1.d-6 
  xp = 0.3d0 
  T = 1.5d0/197.3d0 
  
  nuclei = get_nuclear_ensemble(12,1000,0.3d0,0.6d0,SIZE(nuclei)) 
   
  ! Calculate nuclei with maximum medium correction
  !eos_outside = get_bulk_state(eos_input(1.d-6*nb,xp,T))
 
  !print *,'Calculating nuclei'
  !!$OMP DO 
  !do i=1,SIZE(nuclei)  
  !  nuclei(i) = nucleus_find_properties(nuclei(i)%Z,nuclei(i)%N,nb*xp,eos_outside)  
  !enddo
  !!$OMP END DO
  
  print *,'Calling NSE'  
  eos = get_nse_eos(eos_input(nb,xp,T),nuclei) 
   
end program nse_test 
