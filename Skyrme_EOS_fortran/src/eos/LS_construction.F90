module LS_construction
  use eos_com_mod
  use eos_skyrme_mod   
  implicit none 
  
  real(8), private, parameter :: HBC  = 197.3d0 
  real(8), private, parameter :: MNUC = 939.565d0/HBC
  real(8), private, parameter :: PI   = 3.14592d0 
  ! This definition of ALPHA differs slightly from LS 
  real(8), private, parameter :: ALPHA = 3.d0/(10.d0*MNUC)*(3.d0*PI*PI)**(2.d0/3.d0)
  real(8), private, parameter :: e_ele = SQRT(1.4299764d0/HBC) 
  real(8), private, parameter :: A0 = 60.d0
  real(8), private :: Ss0 = 45.8d0/HBC
  real(8), private :: sigmas = 1.15d0/HBC

contains
  
  function get_LS_state(eosin,include_el) result (eosout)
    use nr
    use nrtype 
    
    type(eos_com), intent(in) :: eosin
    type(eos_com) :: eosout
    logical, intent(in)  :: include_el
    real(8) :: nt,xp,T
    real(8) :: pp,ee,mueq 
    real(8) :: muno,mupo 
    real(8) :: Eel,Pel,mue,usave,nh,xh,nn,xx
    real(8) :: nno,xxo
    real(8) :: uext,fpress 
    real(SP), dimension(5) :: xg,fg
    logical :: bad 
    integer :: ii,nsign 
    real(8) :: ulo,uhi,umid,flo,fhi,fmid,umin,umax
    real(8), dimension(100) :: signlo,signhi 
    
    nt = eosin%nb 
    xp = eosin%xp 
    T  = eosin%T 

    ! Calculate the nuclear part of the EoS      
    !call get_bulk_state(nt,xp,ee,pp,muno,mupo)
    !call get_zero_bulk_mu(nh,xh) 
    xg(1) = LOG(0.14d0*0.5d0)  
    xg(2) = LOG(0.14d0*0.5d0)  
    xg(3) = LOG(1.d-3)
    xg(4) = LOG(1.d-7)
    xg(5) = LOG(1.d-7)
    call newt(xg(1:5),bad,get_f)  
    fg = get_f(xg(1:5)) 
    nh = EXP(xg(1)) + EXP(xg(2))
    print *,'done',bad,EXP(xg(4)) + EXP(xg(5)),nh,EXP(xg(3)),EXP(xg(1))/nh
    !print *,'done',bad,nt,xp,usave,(nt-usave*nh)/(1.d0-usave),nh,fg(1:5)
    
    ! Electron quantities 
    if (include_el) then
      Eel = (3.d0*PI**2*nt*xp)**(4.d0/3.d0)/(4.d0*PI*PI) 
      Pel = Eel/3.d0
      mue = (3.d0*PI*PI*xp*nt)**(1.d0/3.d0) 
    else 
      Eel = 0.d0 
      Pel = 0.d0
      mue = 0.d0
    endif
    
    eosout%ee = (ee + EEl/nt)
    eosout%pp = (pp + Pel)
    eosout%mun = muno
    eosout%mup = mupo
    eosout%mue = mue

  contains
    
    function get_f(x,output) result(ff) 
      real(SP), dimension(:), intent(in) :: x
      logical, optional :: output 
      real(SP), dimension(SIZE(x)) :: ff  
      real(8) :: npi,nni,npo,nno,u
      real(8) :: nto,nti,xpi,xpo,eei,ppi,eeo,ppo,eep
      real(8) :: muni,mupi 
      real(8) :: D,Dp,sig,sigp,beta 
      real(8) :: muscale,pscale,muh,nq,nt_
      type(eos_com)  :: eosi,eoso  
      
      ! Calculate the interior and exterior densities 
      npi = EXP(x(1)) 
      nni = EXP(x(2)) 
      u   = 1.d-4
      nt_ = EXP(x(3))
      if (SIZE(x)>3) then
        npo = EXP(x(4)) 
        nno = EXP(x(5)) 
      else
        npo = MAX(1.d-14,(xp*nt_ - u*npo)/(1.d0-u))
        nno = MAX(1.d-14,((1.d0-xp)*nt_ - u*nno)/(1.d0-u))
      endif 
       
      ! Calculate the interior and exterior bulk properties
      eosi = eos_input_from_densities(npi,nni,T)  
      eoso = eos_input_from_densities(npo,nno,T)  
      eosi = get_bulk_state(eosi)
      eoso = get_bulk_state(eoso)
      
      ! The equilibrium equations
      !call sigma(eosIni%xp,sig,sigp)  
      !call Dsurf(u,D,Dp) 
      beta = 4.5d0*(8.d0*PI/15.d0)**(1.d0/3.d0)*(e_ele*npi*sig)**(2.d0/3.d0) 
      usave = u 
      
      nQ = (Mnuc*T/(2.d0*PI))**1.5d0 
      muh = T*LOG(u*(1.d0-u)*(npi+nni)/(nQ*A0**2.5d0)) 
      muh = 0.d0
       
      ff(1) = (eosi%mun + muh/A0*(1.d0-u) - eoso%mun)/T
      ff(2) = (eosi%mup + muh/A0*(1.d0-u) - eoso%mup)/T 
      ff(3) = (ppi + u*eosi%nb/A0*muh - ppo)/(T*nt_)
      if (SIZE(x)>3) then
        ff(4) = (u*npi + (1.d0-u)*npo)/nt_ - xp 
        ff(5) = (u*nni + (1.d0-u)*nno)/nt_ - (1.d0-xp)
      endif 
       
      if (PRESENT(output)) then 
        if (output) write(6,'(45es11.3)')nno,npo,nni,npi,nt_,u,ff
      endif
      
      if (ANY(ff.ne.ff)) then 
        print *,'f',ff 
        print *,'n',npo,nno,npi,nni,u
        print *,'derived',D,beta,sig,sigp
        print *,'I quit'
        stop  
      endif 
    
    end function get_f
        
  end function get_LS_state 
  
  function get_u_from_pressure_difference(ppi,ppo,bet) result(u) 
    use nr, ONLY: newt 
    use nrtype 
    real(8), intent(in) :: ppi,ppo,bet 
    real(SP) :: u,uu,ul
    real(SP), dimension(1) :: ff,dfdu 
    real(8) :: delta 
    logical :: bad
    real(SP), dimension(1) :: xg 
    integer :: i  
    
    delta = (ppo - ppi)/bet
    
    ! These are the limiting cases 
    !if (delta>5.d0/3.d0) then ! Limits for the actual LS form of D
    if (delta>60.d0) then
      u = 1.d0
      return 
    !elseif (delta<-1.d0/3.d0) then 
    elseif (delta<-1.d0) then 
      u = 0.d0 
      return 
    endif 
    
    ! Do the search if we are inside the limits
    u = 0.5d0 
    do i=1,20 
      uu = u*(1.d0+1.d-6) 
      ul = u*(1.d0-1.d-6) 
      dfdu = (get_f((/uu/)) - get_f((/ul/)))/(uu-ul) 
      ff = get_f((/u/)) 
      if (ff(1)/dfdu(1).ne.ff(1)/dfdu(1)) print *,'looking for u',u,delta,ff(1),dfdu(1) 
      u = MIN(MAX(1.d-2*u,u - ff(1)/dfdu(1)),0.99d0*u)       
    enddo   
     
    !xg(1) = 0.5d0 
    !call newt(xg,bad,get_f)  
    !if (bad) print *,'Find u failed for some reason' 
    !u = xg(1) 
    return 
               
  contains 
  
    function get_f(xx,output) result(ff) 
      real(SP), dimension(:), intent(in) :: xx
      LOGICAL, optional :: output
      real(SP), dimension(SIZE(xx)) :: ff 
      real(8) :: D, Dp  
      call Dsurf(dble(xx(1)),D,Dp)   
      ff(1) = 2.d0/3.d0*D - Dp + (ppi - ppo)/bet
    end function 
  
  end function get_u_from_pressure_difference
   
  subroutine Dsurf(u,D,Dp) ! This is Script{D}/u and Script{D}'  
    real(8), intent(in) :: u 
    real(8), intent(out) :: D,Dp
    real(8) :: a,b,denom,ddeno,uda,mudb 
    a  = 1.d0 - 1.5d0*u**(1.d0/3.d0) + 0.5d0*u + 1.d-50
    uda = -0.5d0*u**(1.d0/3.d0) + 0.5d0*u
    b  = MAX(1.d0 - 1.5d0*(1.d0-u)**(1.d0/3.d0) + 0.5d0*(1.d0-u) + 1.d-50,0.d0)
    mudb = 0.5d0*(1.d0-u)**(1.d0/3.d0) - 0.5d0*(1.d0-u)
    denom = u**2 + (1.d0-u)**2 + 0.6d0*u**2*(1.d0-u)**2 
    ddeno = 2.d0*u - 2.d0*(1.d0 - u) + 1.2d0*u*(1.d0-u)**2 - 1.2d0*u**2*(1.d0-u)  
    
    D  = (1.d0-u)*((1.d0-u)*a**(1.d0/3.d0) + u*b**(1.d0/3.d0))/denom
    
    Dp = (1.d0-2.d0*u)*((1.d0-u)*a**(1.d0/3.d0) + u*b**(1.d0/3.d0))/denom 
    Dp = Dp + ((1.d0-u)**2*a**(-2.d0/3.d0)/3.d0*uda + u**2*b**(-2.d0/3.d0)/3.d0*mudb)/denom
    Dp = Dp + u*(1.d0-u)*(-a**(1.d0/3.d0) + b**(1.d0/3.d0))/denom 
    Dp = Dp - u*(1.d0-u)*((1.d0-u)*a**(1.d0/3.d0) + u*b**(1.d0/3.d0))/denom**2*ddeno
    
    if (D.ne.D) print *,'Bad D',D,u,a,b,denom 
    if (Dp.ne.Dp) print *,'Bad Dp',Dp,u,a,b,uda,mudb,denom,ddeno
  end subroutine Dsurf 

  subroutine sigma(xp,sig,sigp)  
    real(8), intent(in)  :: xp 
    real(8), intent(out) :: sig,sigp
    real(8) :: q,r0
    r0 = (3.d0/(4.d0*PI*0.155d0))**(1.d0/3.d0) 
    q = 384.d0*PI*r0**2*sigmas/Ss0 - 16.d0  
    sig  = sigmas*(16.d0 + q)/(q + xp**(-3) + (1.d0-xp)**(-3)) 
    sigp = sigmas*(16.d0 + q)/(q + xp**(-3) + (1.d0-xp)**(-3))**2*(-3.d0/xp**4 + 3.d0/(1.d0-xp)**4.d0) 
    if (xp>1.d0-1.d-10) sigp = sigmas*(16.d0 + q)*(1.d0-xp)**2*3.d0
    if (xp<1.d-10) sigp = sigmas*(16.d0 + q)*xp**2*(-3.d0)
    if (sig.ne.sig) print *,'sig',q,xp 
    if (sigp.ne.sigp) print *,'sigp',q,xp 
  end subroutine sigma 

end module LS_construction

