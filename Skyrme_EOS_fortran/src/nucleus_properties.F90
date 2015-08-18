program nucleus_properties
  use eos_nse_mod
  use eos_com_mod 
  use eos_skyrme_mod
  use eos_gibbs_mod
  implicit none

  real(8) :: Z,N,ne,vnuc,T,nb,ye
  type(eos_com) :: eos_out, eos_gibbs, eos_low, eos_hi
  type(nucleus) :: nuc
  type(gibbs_region) :: gibbs  
  real(8), parameter :: HBC = 197.3d0 
  integer :: i,j
  
  T = 1.d0/HBC
  nb = 1.d-2 
  ye = 0.4d0 
  
  ! Do the gibbs construction for this region 
  call gibbs%boundary(T)  
  eos_gibbs = gibbs%get_point(eos_input(nb,ye,T)) 
  eos_low = get_bulk_state(eos_input_from_densities(eos_gibbs%npl,eos_gibbs%nnl,T)) 
  eos_hi  = get_bulk_state(eos_input_from_densities(eos_gibbs%npu,eos_gibbs%nnu,T)) 
   
  ne = nb*ye
  do i=100,500
    do j=28,28  
      Z = DBLE(i)
      N = DBLE(j) 
      N = DBLE(i)*(1.d0-eos_low%xp) 
      Z = DBLE(i)*eos_low%xp
      if (i==100) nuc%v = (N+Z)/eos_low%nb
      nuc = nucleus_find_properties(Z,N,ne,eos_hi,nuc%v)  
      write(6,'(20es12.3)')nuc%Z,nuc%N,nuc%v,(nuc%Z+nuc%N)/nuc%v,HBC*nuc%Etot/nuc%A 
    enddo 
  enddo

end program
