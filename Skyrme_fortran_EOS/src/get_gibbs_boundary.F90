program get_gibbs_boundary
  use eos_gibbs_mod
  use eos_com_mod 
  implicit none 

  type(gibbs_region) :: phase_bound 
  type(eos_com) :: eos_out
  real(8) :: T,nb,ye
  integer :: i 
  character(80) :: fname 
  
  T = 2.d0/197.3d0
  call phase_bound%boundary(T) 
  
  ye = 0.3d0 
  do i=1,1001 
    nb = 10.d0**(DBLE(i-1)/1000.d0*(LOG10(1.7d-1) - LOG10(1.d-6)) + LOG10(1.d-6))  
    eos_out = phase_bound%get_point(eos_input(nb,ye,0.4d0*T))
    write(6,'(15es12.3)')eos_out%nb,eos_out%pp,eos_out%nnl,eos_out%npl, &
        eos_out%nnu,eos_out%npu,eos_out%u 
  enddo 
  
   
  write(fname,'(a,f0.1,a)')'PhaseBoundary.T',T*197.3d0,'MeV'
  open(12,file=TRIM(fname))  
  do i=1,phase_bound%npoints
      write(12,'(50es16.3e3)')phase_bound%mun(i),phase_bound%npl(i), &
                              phase_bound%nnl(i),phase_bound%npu(i), &
                              phase_bound%nnu(i),phase_bound%pp(i),  &
                              phase_bound%mup(i) 
  enddo    
  close(12) 
   
   
    
end program get_gibbs_boundary 
