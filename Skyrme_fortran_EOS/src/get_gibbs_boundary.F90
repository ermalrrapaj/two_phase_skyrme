program get_gibbs_boundary
  use eos_gibbs_mod, ONLY: phase_bound
  use eos_com_mod, ONLY: eos_com 
  implicit none 

  type(gibbs_region) :: phase_bound 
  type(eos_com) :: eos_out
  real(8) :: T,nb,ye
  integer :: i 
  character(80) :: fname 
 
   
  ! Use the gibbs_region module to find the phase boundary at fixed temperature
  T = 2.d0/197.3d0
  call phase_bound%boundary(T) 
  
  ye = 0.3d0 
  do i=1,1001 
    nb = 10.d0**(DBLE(i-1)/1000.d0*(LOG10(1.7d-1) - LOG10(1.d-6)) + LOG10(1.d-6))  
    eos_out = phase_bound%get_point(eos_input(nb,ye,0.4d0*T))
    write(6,'(15es12.3)')eos_out%nb,eos_out%pp,eos_out%nnl,eos_out%npl, &
        eos_out%nnu,eos_out%npu,eos_out%u 
  enddo 
  
  ! Write out the phase boundary data to file 
  write(fname,'(a,f0.1,a)')'PhaseBoundary.T',T*197.3d0,'MeV'
  open(12,file=TRIM(fname))  
  write(12,'(a)') '# [1] mu_n (fm^-1)'
  write(12,'(a)') '# [2] n_{p,low density} (fm^-3)'
  write(12,'(a)') '# [3] n_{n,low density} (fm^-3)'
  write(12,'(a)') '# [4] n_{p,high density} (fm^-3)'
  write(12,'(a)') '# [5] n_{n,high density} (fm^-3)'
  write(12,'(a)') '# [6] Pressure (fm^-4)'
  write(12,'(a)') '# [7] mu_p (fm^-1)'
  do i=1,phase_bound%npoints
      write(12,'(50es16.3e3)')phase_bound%mun(i),phase_bound%npl(i), &
                              phase_bound%nnl(i),phase_bound%npu(i), &
                              phase_bound%nnu(i),phase_bound%pp(i),  &
                              phase_bound%mup(i) 
  enddo    
  close(12) 
   
   
    
end program get_gibbs_boundary 
