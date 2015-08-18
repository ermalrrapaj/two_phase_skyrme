module eos_skyrme_mod
  use eos_com_mod
   
  implicit none 
  
  ! Constants 
  real(8), private, parameter :: HBC  = 197.3d0 
  real(8), private, parameter :: MNUC = 939.565d0/HBC
  real(8), private, parameter :: PI   = 3.14592d0 
  ! This definition of ALPHA differs slightly from LS 
  real(8), private, parameter :: ALPHA = 3.d0/(10.d0*MNUC)*(3.d0*PI*PI)**(2.d0/3.d0)
  real(8), private, parameter :: e_ele = SQRT(1.4299764d0/HBC) 
  
  ! Constants 
  real(8), private :: ns0 = 0.155d0 
  real(8), private :: Sv0 = 29.3d0/HBC
  real(8), private :: Ss0 = 45.8d0/HBC
  real(8), private :: sigmas = 1.15d0/HBC
  
  ! Parameters of the force
  real(8), private :: a = -285.1d0/HBC 
  real(8), private :: b = -107.1d0/HBC 
  real(8), private :: c = 968.d0/HBC 
  real(8), private :: d = 0.d0

  real(8), private :: f = 0.d0
  real(8), private :: g = 0.d0
  
  real(8), private :: delta = 2.002d0 
  
contains
  
  subroutine set_Skyrme(ns,BE,K,Sv,mstarom,L,use_Lin,Ks,use_Ksin) 
    ! Set the skyrme parameters using properties of matter at saturation 
    use nr 
    use nrtype 
   
    real(8), intent(in) :: ns,BE,K,Sv,mstarom,L,Ks
    logical, intent(in) :: use_Lin,use_Ksin
    real(8) :: alphaLS,ft
    real(SP), dimension(6) :: xg,fg
    logical(LGT) :: bad
    logical :: use_L,use_Ks
    
    ns0 = ns 
    Sv0 = Sv/HBC 
     
    ! Set the effective mass parameter and the kinetic energy parameter
    ft = (1.d0/mstarom - 1.d0)
    alphaLS = ALPHA*(ns/2.d0)**(2.d0/3.d0)/mstarom  

    ! Set the initial guess 
    xg(1) = a*ns 
    xg(2) = b*ns
    xg(3) = c*ns**delta   
    xg(4) = delta
    xg(5) = d*ns**delta  
    xg(6) = g*ns
     
    ! Solve for parameters
    use_L = use_Lin 
    use_Ks = use_Ksin  
    call newt(xg,bad,get_f)  
    if (bad) then
      print *,'Skyrme parameter search failed' 
      stop
    endif
    
    ! Transform back to LS style parameters   
    a = xg(1)/ns 
    b = xg(2)/ns
    delta = xg(4) 
    c = xg(3)/ns**delta 
    f = ft/ns 
    d = xg(5)/ns**delta  
    g = xg(6)/ns 
  
  contains 
    
    function get_f(x,output) result(ff) 
      real(SP), dimension(:), intent(in) :: x 
      logical, optional :: output
      real(SP), dimension(SIZE(x)) :: ff 
      real(8) :: ta,tb,tc,td,tg,tdel
      ta = 0.d0; tb = 0.d0; tc = 0.d0; tdel = 0.d0; td = 0.d0; tg = 0.d0 
      
      if (SIZE(x)>0) ta = x(1)  
      if (SIZE(x)>1) tb = x(2)  
      if (SIZE(x)>2) tc = x(3)
      if (SIZE(x)>3) tdel = x(4)  
      if (SIZE(x)>4) td = x(5)  
      if (SIZE(x)>5) tg = x(6)  
      
      if (SIZE(x)>0) ff(1) = (10.d0/9.d0*alphaLS*(1.d0 + 4.d0*ft) + (2.d0*(ta + tb) + (tc + td)*tdel*(1 + tdel)))*HBC*9.d0/K-1.d0
      if (SIZE(x)>1) ff(2) = (alphaLS*(2.d0 + 5.d0*ft)/3.d0 + ta + tb + (tc + td)*tdel)/alphaLS 
      if (SIZE(x)>2) ff(3) = (alphaLS*(1.d0+ft) + ta + tb + tc + td)*HBC/BE + 1.d0
      if (SIZE(x)>3) ff(4) = (alphaLS*((2.d0**(2.d0/3.d0)-1.d0)*(1.d0 + ft) - 1.d0*2.d0**(2.d0/3.d0)*tg) - (tb+td))*HBC/Sv - 1.d0
      if (use_L) then
        if (SIZE(x)>4) ff(5) = (alphaLS*((2.d0**(2.d0/3.d0)-1.d0)*(2.d0 + 5.d0*ft) - 5.d0*2.d0**(2.d0/3.d0)*tg) & 
                               - 3.d0*(tb + td*tdel)) - L/HBC
      else
        if (SIZE(x)>4) ff(5) = td 
      endif
      if (use_Ks .and. use_L) then  
        if (SIZE(x)>5) ff(6) = -Ks/HBC + (9*td*(-1 + tdel)*tdel - 2*alphaLS*(1 - 2**(2.d0/3.d0) &
                               +5*(-1 + 2**(2.d0/3.d0))*ft - 5*2**(2.d0/3.d0)*tg))/3.d0
      else 
        if (SIZE(x)>5) ff(6) = tg
      endif 
    end function get_f 

  end subroutine set_Skyrme 

  subroutine set_Skyrme_LS(ns,BE,K,Sv) 
   ! Set the Skyrme parameters only using the binding energy, compressibility, 
   ! and symmetry energy
   real(8), intent(in) :: ns,BE,K,Sv
   real(8) :: alphaLS 

   alphaLS = ALPHA*(ns/2.d0)**(2.d0/3.d0)  
   
   delta = (K/HBC + 2.d0*alphaLS)/(9.d0*BE/HBC + 3.d0*alphaLS) 
   b = (alphaLS*(2.d0**(2.d0/3.d0)-1.d0) - Sv/HBC)/ns   
   a = (delta*(alphaLS+BE/HBC) - 2.d0/3.d0*alphaLS)/(ns*(1.d0-delta)) - b
   c = (K/HBC + 2.d0*alphaLS)/(9.d0*delta*(delta-1.d0)*ns**delta) 
  end subroutine  
  
  subroutine get_zero_bulk_mu(nt,xp) 
    ! Find were the neutron and proton chemical potentials are both zero
    use nrtype
    use nr 
    real(8), intent(out) :: xp, nt
    real(SP), dimension(2) :: xg  
    logical :: check 
    xg(1) = 0.155d0 
    xg(2) = 0.5d0
    call newt(xg,check,get_f) 
    nt = DBLE(xg(1)) 
    xp = DBLE(xg(2)) 
  contains
    function get_f(xx,output) result(ff) 
      real(SP), intent(IN), dimension(:) :: xx 
      logical, optional :: output
      real(SP), dimension(size(xx)) :: ff
      real(8) :: ee,pp,mun,mup  
      type(eos_com) :: eosin 
      type(eos_com) :: eosout 
      eosin = eos_input(nb=xx(1),xp=xx(2))
      eosout = get_bulk_state(eosin)
      ff(1) = eosout%mun
      ff(2) = eosout%mup 
    end function get_f 
  end subroutine get_zero_bulk_mu 
  
  function get_bulk_state_from_mixed(eos_in) result(eos_out)
    ! Given the proton number density and the neutron chemical potential,
    ! find the neutron number density and return the eos properties
    use nrtype 
    use nr 
    use bisect_mod 
    type(eos_com), intent(in) :: eos_in 
    type(eos_com) :: eos_out
    logical :: bad
    real(8) :: xx,yy,yu,yl,dydx,nbl,nbh
    integer :: i  
    
    !xx = LOG(eos_in%np)
    !call newt(xx,bad,get_f) 
    nbl = 1.d-8
    nbh = 1.d0
    yl = get_f(LOG(nbl))  
    yu = get_f(LOG(nbh))  
    do i=1,50
      ! Check that we aren't too far away
      if (yl<-HUGE(yl)) then
        nbl = nbl*1.99d0  
        yl = get_f(LOG(nbl))
      endif 

      if( yu>HUGE(yu)) then
        nbh = nbh/1.99d0  
        yu = get_f(LOG(nbh)) 
      endif 
      
      if (yl*yu<=0.d0 .and. yu<HUGE(yu) .and. yl>-HUGE(yl)) then
        exit
      endif 
      
      ! We are the same sign if we made it here
      if (yl>0.d0) then 
        nbl = nbl/2.d0
        yl = get_f(LOG(nbl))
      endif
       
      if (yu<0.d0) then
        nbh = nbh*2.d0
        yu = get_f(LOG(nbh)) 
      endif   
       
    enddo 
    
    if (i>50) then
      !print *,'Could not find interval with zero' 
      !print *,'yl',nbl,yl,nbh,yu
      return
    endif  
     
    xx = bisect(LOG(nbl),LOG(nbh),bad,get_f) 
    if (.not. bad) then
      eos_out = get_bulk_state(eos_input_from_densities(eos_in%np,EXP(xx),eos_in%T))
    endif  
  contains 
   
    function get_f(x_) result(y_)
      real(SP), intent(IN) :: x_
      real(SP) :: y_
      type(eos_com) :: eoso_ 
      eoso_ = get_bulk_state(eos_input_from_densities(eos_in%np,EXP(x_),eos_in%T))  
      y_ = (eoso_%mun - eos_in%mun)/MIN(ABS(eos_in%mun),eos_in%T)
    end function get_f     
    
  end function get_bulk_state_from_mixed 
   
  function get_bulk_state_from_mu(eos_in) result (eos_out)  
    ! Given the neutron and proton chemical potentials, find the neutron and 
    ! proton number densities and return the eos properties
    use nrtype 
    use nr 
    type(eos_com), intent(in) :: eos_in 
    type(eos_com) :: eos_out(2)
    logical :: bad
    real(8), SAVE :: xlo(2) = LOG(1.d-8)
    real(8), SAVE :: xhi(2) = LOG(1.d1)
    
    eos_out(1)%full_set = .false. 
    eos_out(2)%full_set = .false. 
    
    if (.not.eos_in%cpotentials_set) then
      print *,'Passing bad input to get bulk state from mu' 
      stop
    endif 

    xlo = LOG(1.d-8) ! This is maybe not the best way to do this
    call newt(xlo,bad,get_f)  
    if (.not. bad) then
      eos_out(1) = get_bulk_state(eos_input_from_densities(EXP(xlo(1)),EXP(xlo(2)),eos_in%T))
    endif  
    
    !xx = LOG(3.d-2)   
    !call newt(xx(1:2),bad,get_f)  
    !if (.not. bad) then
    !  eos_out = get_bulk_state(eos_input_from_densities(EXP(xx(1)),EXP(xx(2)),eos_in%T))
    !  write(6,'(10es11.2)')eos_out%np,eos_out%nn,eos_out%mup/eos_out%T,eos_out%mun/eos_out%T,eos_out%pp 
    !endif  
    
    xhi = LOG(1.d1)   
    call newt(xhi,bad,get_f)  
    if (.not. bad) then
      eos_out(2) = get_bulk_state(eos_input_from_densities(EXP(xhi(1)),EXP(xhi(2)),eos_in%T))
    endif
    
    ! Check that the return isn't double valued and the actual one is in the first space 
    if (eos_out(2)%full_set .and. eos_out(1)%full_set) then
      if (ABS(eos_out(1)%np/(eos_out(2)%np)-1.d0) < 1.d-6 .and. ABS(eos_out(1)%nn/(eos_out(2)%nn)-1.d0) < 1.d-6) then
        eos_out(2)%full_set = .false. 
      endif
    elseif (eos_out(2)%full_set) then
      eos_out(1) = eos_out(2) 
      eos_out(2)%full_set = .false.    
    endif 
       
  contains 
    
    function get_f(x_,output_) result(y_)
      real(SP), intent(IN), dimension(:) :: x_
      logical, optional :: output_
      real(SP), dimension(size(x_)) :: y_
      type(eos_com) :: eoso_ 
      eoso_ = get_bulk_state(eos_input_from_densities(EXP(x_(1)),EXP(x_(2)),eos_in%T))
      y_(1) = (eoso_%mun - eos_in%mun)/MAX(ABS(eos_in%mun),1.d-10) 
      y_(2) = (eoso_%mup - eos_in%mup)/MAX(ABS(eos_in%mup),1.d-10)
    end function get_f     
  
  end function get_bulk_state_from_mu 

  function get_bulk_state(eos_in) result (eos_out)  
    ! For given baryon number density and proton fraction, find the 
    ! thermodynamic properties
    type(eos_com), intent(in) :: eos_in 
    type(eos_com) :: eos_out

    real(8) :: nt,xp
    real(8) :: pp,ee,mun,mup
    
    real(8) :: ET,EU  
    real(8) :: PT,PU  
    real(8) :: xn53,xp53,nt53,np,nn,taun,taup
    real(8) :: momsp, momsn, Vp, Vn, mue, etap, etan,T
    double precision :: ifermi32,zfermi32
    double precision :: ifermi12,zfermi12
    
    nt = MIN(eos_in%nb,1.d4)
    xp = MIN(MAX(eos_in%xp,0.d0),1.d0) 
    T  = eos_in%T
     
    ! Useful quantities 
    xn53 = (1.d0-xp)**(5.d0/3.d0) 
    xp53 = xp**(5.d0/3.d0) 
    nt53 = nt**(5.d0/3.d0)
    nn = (1.d0 - xp)*nt 
    np = xp*nt
    
    ! Ratio of vacuum mass over effective mass
    momsp = (1.d0 + f*(nn+np) + g*(nn-np))
    momsn = (1.d0 + f*(nn+np) - g*(nn-np))
    
    ! Single particle potentials
    etap = ifermi12(2.d0*PI**2*np*(2.d0*Mnuc/momsp*T)**(-1.5d0)) 
    etan = ifermi12(2.d0*PI**2*nn*(2.d0*Mnuc/momsn*T)**(-1.5d0)) 
    
    taup = (2.d0*Mnuc/momsp*T)**(5.d0/2.d0)*zfermi32(etap)/(2.d0*PI**2) 
    taun = (2.d0*Mnuc/momsn*T)**(5.d0/2.d0)*zfermi32(etan)/(2.d0*PI**2) 
    
    Vp = (taup*(f - g) + taun*(f + g))/(2.d0*Mnuc) & 
        + 2.d0*a*nt + 4.d0*b*nn + c*(1.d0+delta)*nt**delta &
        + 4.d0*d*nt**(delta-2.d0)*(nn*nt + (delta-1.d0)*nn*np)  
    Vn = (taup*(f + g) + taun*(f - g))/(2.d0*Mnuc) & 
        + 2.d0*a*nt + 4.d0*b*np + c*(1.d0+delta)*nt**delta &
        + 4.d0*d*nt**(delta-2.d0)*(np*nt + (delta-1.d0)*nn*np)  
     
    mup = etap*T + Vp
    mun = etan*T + Vn
     
    ! Thermodynamic quantities
    ET = taup*momsp/(2.d0*Mnuc) + taun*momsn/(2.d0*Mnuc) 

    PT = ((5.d0/6.d0*momsp - 0.5d0)*taup + (5.d0/6.d0*momsn - 0.5d0)*taun)/Mnuc
    
    EU = (a + 4.d0*b*xp*(1.d0-xp))*nt**2 + (c + 4.d0*d*xp*(1.d0-xp))*nt**(1.d0+delta)
     
    PU = (a + 4.d0*b*xp*(1.d0-xp))*nt**2 + (c + 4.d0*d*xp*(1.d0-xp))*delta*nt**(1.d0+delta) 
    
    eos_out = eos_in
    eos_out%ss = (5.d0*taup/(6.d0*Mnuc*T)*momsp - np*etap + 5.d0*taun/(6.d0*Mnuc*T)*momsn - nn*etan)/nt
    eos_out%pp = PT + PU
    eos_out%ee = (ET + EU)/nt
    eos_out%mun = mun  
    eos_out%mup = mup
    eos_out%mue = 0.d0 
    eos_out%full_set = .true.
    !pp = np*mup + nn*mun - ee*nt
      
    if (eos_out%pp.ne.eos_out%pp ) then
      print *,'bulk nt: ',nt,' xp: ',xp,' T: ',T,' PT:',PT,' PU:',PU
      stop
    endif    
    
    return 
  end function get_bulk_state   

end module eos_skyrme_mod 

