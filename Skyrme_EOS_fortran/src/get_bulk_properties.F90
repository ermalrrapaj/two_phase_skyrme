program test 
  
  use eos_skyrme_mod
  use eos_com_mod
  use bisect_mod
  use nr
  use nrtype
  
  implicit none 

  real(8) :: nt,xp,ee,pp,mueq,rho,e1,e2,mun,mup,musym,pps,ees
  real(8) :: eoff = 8.85e18
  real*8,parameter :: rho_gf   = 1.61930347d-18
  real*8,parameter :: press_gf = 1.80171810d-39
  real*8,parameter :: eps_gf   = 1.11265006d-21
  
  real(8) :: u,f,nlo,nhi,nloo,nhio,npl,eta,T,nhi0,ptransition
  real(SP) :: xx(4),xo(4),xg(4)
  integer :: np,i,j
  logical :: checklo,checkhi
  double precision :: zfermi32,zfermi52 
  type(eos_com) :: eosout 
  
  !call set_Skyrme(0.155d0,16.d0,375.d0,29.3d0,0.95d0,55.0d0,.true.,37.d0,.true.) 
  call set_Skyrme(0.155d0,16.d0,375.d0,29.3d0,1.d0,55.0d0,.false.,37.d0,.false.) 
  
  ! Find the forbidden region 
  T = 5.0d0/197.3d0 
  do i=1,1
    xp   = 10.d0**(LOG10(0.5d0) +  (LOG10(1.d-3)-LOG10(0.5d0))*DBLE(i-1)/999.d0)
    if (i>1) then
      nlo = bisect(nloo,nhio*0.99d0,checklo,zero_pressure) 
      nhi = bisect(nlo*1.0001d0,nhio,checkhi,zero_pressure) 
    else 
      nlo = bisect(1.0d-5,0.1d0,checklo,zero_pressure) 
      nhi = bisect(0.1d0,0.8d0,checkhi,zero_pressure) 
      xg(1) = LOG(xp*nlo*0.7d0)
      xg(2) = LOG((1.d0-xp)*nlo*0.8d0)
      xg(3) = LOG(xp*nhi*0.9d0)
      xg(4) = LOG((1.d0-xp)*nhi)
      nhi0 = nhi
    endif
    if (checklo .or. checkhi) exit
    !print *,xp,nlo,nhi,nlo*xp,nlo*(1.0-xp),nhi*xp,nhi*(1.0-xp)
    nloo = nlo 
    nhio = nhi
  enddo
  
  xx = xg 
  ptransition = 0.1d0/197.3d0
  call newt(xx(1:4),checklo,phase_boundary)  
  ptransition = ptransition*0.33d0
  call newt(xx(1:4),checklo,phase_boundary)  
  ptransition = ptransition*0.33d0
  call newt(xx(1:4),checklo,phase_boundary)  
  ptransition = ptransition*0.33d0
  call newt(xx(1:4),checklo,phase_boundary)  
  ptransition = ptransition*0.6d0
  call newt(xx(1:4),checklo,phase_boundary)  
  ptransition = ptransition*0.8d0
  call newt(xx(1:4),checklo,phase_boundary)  
  T = 4.5d0/197.3d0
  call newt(xx(1:4),checklo,phase_boundary)  
  T = 1.d0/197.3d0
  call newt(xx(1:4),checklo,phase_boundary)  
  T = 0.5d0/197.3d0
  call newt(xx(1:4),checklo,phase_boundary)  
  
  do i=1,500
    xo = xx
    call newt(xx(1:4),checklo,phase_boundary)  
    if (checklo) then
      ptransition = ptransition*(1.d0+1.d-2)
      xx = xo 
      cycle
    endif 
     
    print *,ptransition*197.3d0,T*197.3d0,EXP(xx(1)),EXP(xx(2)),EXP(xx(1))/(EXP(xx(1))+EXP(xx(2))),& 
          EXP(xx(3)),EXP(xx(4)),EXP(xx(3))/(EXP(xx(3))+EXP(xx(4))),checklo
    ptransition = ptransition/(1.d0+2.d-2)
    !T = T/(1.d0+2.d-2)
  enddo 
  return 
   
  xx = xg 
  ptransition = 0.1d0/197.3d0
  do i=1,100
    call newt(xx(1:4),checklo,phase_boundary)  
    if (checklo) exit 
    print *,ptransition*197.3d0,EXP(xx(1)),EXP(xx(2)),EXP(xx(1))/(EXP(xx(1))+EXP(xx(2))),& 
          EXP(xx(3)),EXP(xx(4)),EXP(xx(3))/(EXP(xx(3))+EXP(xx(4))),checklo
    ptransition = ptransition*(1.d0+1.d-2)
  enddo 
  return
  !xx(3) = xx(2)
  !return
  !npl = EXP(xx(1)) 
  !call newt(xx(1:3),checklo,phase_boundary)  
  
  do i=1,1000  
    nt = 10.d0**(LOG10(0.2d0) + (LOG10(1.d-4) - LOG10(0.2d0))*DBLE(i)/999.d0) 
    xp = 0.5 
    T  = 1.d0/197.3d0
    eosout = get_bulk_state(eos_input(nt,xp,T))
    print *,nt,xp,T,eosout%pp,eosout%mun,eosout%ee 
  enddo
     
  return
   
  !do i=1,1000
  !  nt = 10.d0**(LOG10(0.2d0) + (LOG10(1.d-4) - LOG10(0.2d0))*DBLE(i)/999.d0) 
  !  xp = 0.5d0
  !  write(6,'(2es12.3)',ADVANCE="no") nt,xp
  !  do j=1,10
  !    xp = DBLE(j-1)/18.d0
  !    call get_bulk_state(nt,xp,ee,pp,mun,mup) 
  !    write(6,'(es12.3)',ADVANCE="no")pp
  !  enddo   
  !  write(6,'(a)')'' 
  !enddo
   
contains
 
  function zero_pressure(n_) result(p_) 
    implicit none 
    real(8) :: p_ 
    real(8), intent(in) :: n_
    real(8) :: e_,mun_,mup_
    type(eos_com) :: eosin 
    type(eos_com) :: eosout
    eosin = eos_input(n_,xp,T)
    eosout = get_bulk_state(eosin) 
    p_ = eosout%pp  
  end function zero_pressure 
  
  function phase_boundary(x_,out_) result(y_)
    implicit none 
    real(SP), dimension(:), intent(in) :: x_
    logical, optional :: out_
    real(SP), dimension(SIZE(x_)) :: y_  
    type(eos_com) :: eosl_,eosh_
    
    ! Set the state of the two phases
    eosl_ = eos_input_from_densities(EXP(DBLE(x_(1))),EXP(DBLE(x_(2))),T) 
    eosh_ = eos_input_from_densities(EXP(DBLE(x_(3))),EXP(DBLE(x_(4))),T) 
    eosl_ = get_bulk_state(eosl_)  
    eosh_ = get_bulk_state(eosh_)  
    
    ! How close are we to equilibrium
    y_(1) = eosl_%pp/ptransition - 1.d0
    y_(2) = eosh_%pp/ptransition - 1.d0
    y_(3) = (eosh_%mun - eosl_%mun)/(eosh_%mun + eosl_%mun)
    y_(4) = (eosh_%mup - eosl_%mup)/(eosh_%mup + eosl_%mup)

  end function phase_boundary

end program
