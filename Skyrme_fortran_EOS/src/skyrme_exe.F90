program test 

  use eos_skyrme_mod
  
  real(8) :: nt,xp,ee,pp,mueq,rho,e1,e2,mun,mup,musym,pps,ees
  real(8) :: eoff = 8.85e18
  real*8,parameter :: rho_gf   = 1.61930347d-18
  real*8,parameter :: press_gf = 1.80171810d-39
  real*8,parameter :: eps_gf   = 1.11265006d-21
  
  real(8) :: u,f
  integer :: np,i 
  
  !call set_Skyrme(0.155d0,16.d0,375.d0,29.3d0,0.95d0,55.0d0,.true.,37.d0,.true.) 
  call set_Skyrme(0.155d0,16.d0,375.d0,29.3d0,1.d0,55.0d0,.false.,37.d0,.false.) 
  do i=1,1000
    nt = 10.d0**(LOG10(0.2d0) + (LOG10(1.d-4) - LOG10(0.2d0))*DBLE(i)/999.d0) 
    xp = 0.5d0
    !xp = 0.01d0*DBLE(i-1)
    !nt = 1.8d-1 
    xp = 0.5d0
    call get_bulk_state(nt,xp,ees,pps,musym,mup) 
    xp = 0.d0
    call get_bulk_state(nt,xp,ee,pp,mun,mup) 
    print *,nt,xp,musym,mun,mup,pp,pps,ee,ees 
    !call get_state(8.d-2/DBLE(i),5.d-1,ee,pp,mueq,.true.) 
  enddo 
  return 

  write(6,'(es14.5)') eoff
   
  np = 500  
  do i=1,np
    rho = 10.d0**(6.d0 + dble(i-1)/dble(np-1)*9.5d0)
    nt = 6.0221367d23*1.d-39*rho 
    xp = get_equilibrium_xp(nt)  
    call get_state(nt,xp,ee,pp,mueq,.true.) 
    call get_state(nt,0.5d0,e1,pp,mueq,.false.)
    call get_state(nt,0.d0,e2,pp,mueq,.false.)
    write(6,'(10es20.8)') log10(rho),log10(pp*1.6021772d-6*1.0d39*press_gf), &
      log10((ee*1.6021772d-6*6.0221367d23+eoff)*eps_gf), xp, e2-e1
    !write(6,'(10es14.5)') nt,e2-e1,xp
  enddo

end program
