module eos_nse_mod
  
  use eos_com_mod 
  use bisect_mod 
  use eos_skyrme_mod
  use nrtype, ONLY: SP 
  use nr 

  implicit none 
  
  real(8), private, parameter :: PI = 3.14159d0 
  real(8), private, parameter :: HBC = 197.3d0  
  real(8), private, parameter :: Mnuc = 939.d0/197.3d0
  
  type :: nucleus 
    real(8) :: N,Z,A
    real(8) :: v,Etot,Ec,Es,F
  end type 
   
contains 

  function get_nuclear_ensemble(Amin,Amax,Yemin,Yemax,nmax) result(nuclei) 
    integer, intent(in) :: Amin,Amax,nmax
    real(8), intent(in) :: Yemin,Yemax
    type(nucleus), dimension(nmax) :: nuclei  
    integer :: i,j,Zmin,Zmax,idx
    idx = 1 
    do i=Amin,Amax
      Zmin = FLOOR(Yemin*DBLE(i))
      Zmax = CEILING(Yemax*DBLE(i))
      do j=Zmin,Zmax
        nuclei(idx)%A = DBLE(i)  
        nuclei(idx)%Z = DBLE(j)
        nuclei(idx)%N = nuclei(idx)%A - nuclei(idx)%Z
        idx = idx + 1 
        if (idx>nmax) exit  
      enddo   
      if (idx>nmax) exit  
    enddo 
    return   
  end function get_nuclear_ensemble 
   
  function get_nse_eos(eos_in,nuclei) result(eos_fin)
    type(eos_com), intent(in)  :: eos_in
    type(nucleus), dimension(:), intent(inout) :: nuclei 
    type(eos_com) :: eos_fin 
    
    logical :: bad 
    integer :: i  
    real(8) :: T,xx(2),yy(2),err,nn,np,u 
    type(eos_com) :: eos_bulk
    logical :: do_update 

    T = eos_in%T 
    do_update = .true. 
     
    xx(1) = LOG(1.d-8*eos_in%np)  
    xx(2) = LOG(1.d-8*eos_in%nn)  
    call newt(xx,bad,nse_zero)  
    yy = nse_zero(xx) 
    err = SUM(ABS(yy)) 
    
    if (err>1.d-8 .or. (err.ne.err)) then
      ! Ok, the first try failed so try stepping down in temperature
      T = MAX(8.d0/HBC,eos_in%T)
      xx(1) = LOG(1.d-8*eos_in%np)  
      xx(2) = LOG(1.d-8*eos_in%nn)
      do i=1,100
        print *,i,HBC*T,EXP(xx)
        call newt(xx,bad,nse_zero)
        if (bad) return 
        if (T==eos_in%T) exit 
        T = MAX(T/1.05d0,eos_in%T)   
      enddo    
      if (i>100) return  
    endif  
    
    do_update = .true.  
    call newt(xx,bad,nse_zero)  
     
    eos_bulk = get_bulk_state(eos_input_from_densities(EXP(xx(1)),EXP(xx(2)),T)) 
    call nucleus_get_densities(nuclei,eos_bulk,eos_in%np,do_update,nn,np,u)
    print *,EXP(xx),nn,np,u
      
  contains 
    function nse_zero(x_,output) result(y_) 
      real(SP), dimension(:), intent(in) :: x_
      logical, optional :: output  
      real(SP), dimension(SIZE(x_)) :: y_
      type(eos_com) :: eos_bulk_ 
      real(8) :: nn_,np_,u_  
      eos_bulk_ = get_bulk_state(eos_input_from_densities(EXP(x_(1)),EXP(x_(2)),T)) 
      call nucleus_get_densities(nuclei,eos_bulk_,eos_in%np,do_update,nn_,np_,u_)
      y_(1) = (nn_ + np_ + (1.d0 - u_)*(eos_bulk_%nn + eos_bulk_%np))/eos_in%nb - 1.d0 
      y_(2) = (nn_ - np_ + (1.d0 - u_)*(eos_bulk_%nn - eos_bulk_%np))/eos_in%nb - (1.d0 - 2.d0*eos_in%xp)  
      write(6,'(20es12.3)')y_,nn_,np_,u_,EXP(x_)
    end function nse_zero
  end function get_nse_eos 
    
  subroutine nucleus_get_densities(nuclei,eos_out,ne,do_update,nn,np,unuc) 
    type(nucleus), dimension(:), intent(inout) :: nuclei
    type(eos_com), intent(in) :: eos_out
    logical, intent(in) :: do_update
    real(8), intent(in) :: ne 
    real(8), intent(out) :: nn,np,unuc 
    !real(8), intent(out), dimension(SIZE(nuclei)) :: nnuc
    real(8) :: eta,nQ,nc
    integer :: i 

    nQ = (Mnuc*eos_out%T/(2.d0*PI))**1.5d0
    nn = 0.d0 
    np = 0.d0 
    unuc = 0.d0
     
    !!$OMP DO REDUCTION(+:nn,np,unuc) PRIVATE(eta,nc)
    do i=1,SIZE(nuclei)
      if (do_update) nuclei(i) = nucleus_find_properties(nuclei(i)%Z,nuclei(i)%N,ne,eos_out,nuclei(i)%v) 
      eta = (nuclei(i)%Z*eos_out%mup + nuclei(i)%N*eos_out%mun - eos_out%pp*nuclei(i)%v - nuclei(i)%F)/eos_out%T
      nc = nQ*nuclei(i)%A**1.5d0*EXP(eta) 
      nn   = nn   + nc*nuclei(i)%N
      np   = np   + nc*nuclei(i)%Z
      unuc = unuc + nc*nuclei(i)%v
    enddo    
    !!$OMP END DO 

  end subroutine nucleus_get_densities

  function nucleus_find_properties(Z,N,ne,eos_outside,vg) result(nuc) 
    real(8), intent(in) :: Z,N,ne ! Proton number, neutron number, average electron density
    type(eos_com), intent(in) :: eos_outside ! Properties of the bulk nucleons 
    real(8), intent(in), optional :: vg 
    type(nucleus) :: nuc
    
    real(8) :: pscale,xg,Es(3),Ec(3),vnuc
    logical :: check,do_output = .false.
    type(eos_com) :: bulk 

    if (PRESENT(vg).and.vg>0.d0) then
      xg = LOG(vg) 
    else 
      xg = LOG((N+Z)/0.16d0)
    endif 
    
    if (eos_outside%pp<1.d-4) then
      pscale = 1.d-4
    else 
      pscale = eos_outside%pp 
    endif  
     
    vnuc = EXP(bisect(xg,check=check,usrfunc=pressure_equation))
    
    if (check) then
      print *,'bad point'
      do_output = .true.  
      xg = pressure_equation(xg) 
      stop
    endif    
    
    Es = get_Esurf(Z,vnuc,eos_outside%np,ne)
    Ec = get_Ecoul(Z,vnuc,eos_outside%np,ne)
    bulk = get_bulk_state(eos_input_from_densities(Z/vnuc,N/vnuc,eos_outside%T)) 
    nuc%N = N 
    nuc%Z = Z 
    nuc%A = Z + N 
    nuc%v = vnuc
    nuc%Etot = bulk%ee*(Z+N) + Es(1) + Ec(1)
    nuc%Ec = Ec(1)    
    nuc%Es = Es(1)    
    nuc%F  = nuc%Etot - eos_outside%T*bulk%ss*(Z+N)  
  contains 
    function pressure_equation(x_) result(y_) 
      real(8), intent(IN) :: x_
      real(8) :: y_
      type(eos_com) :: bulk_
      real(8) :: Es_(3),Ec_(3),v_
      v_ = EXP(x_) 
      bulk_ = get_bulk_state(eos_input_from_densities(Z/v_,N/v_,eos_outside%T)) 
      Es_ = get_Esurf(Z,v_,eos_outside%np,ne)
      Ec_ = get_Ecoul(Z,v_,eos_outside%np,ne)
      y_ = (bulk_%pp - Ec_(2) - Es_(2))/pscale - eos_outside%pp/pscale  
      if (do_output) write(6,'(20es12.3)')v_,bulk_%pp,-Ec_(2),-Es_(2),eos_outside%pp,y_,eos_outside%np,bulk_%np,ne
    end function pressure_equation 
  end function nucleus_find_properties
  
  function get_Esurf(Z,v,npo,ne) result(Es) 
    real(8), intent(in) :: Z,v,npo,ne 
    real(8), dimension(3) :: Es
    real(8), parameter :: sigmas = 1.15d0/197.3d0 
    real(8) :: rnuc 

    Es(3) = 0.d0 
    rnuc = (3.d0*v/(4.d0*PI))**(1.d0/3.d0) 
    Es(1) = 4.d0*PI*sigmas*rnuc**2 
    Es(2) = 2.d0/3.d0*Es(1)/v 
     
  end function get_Esurf 
   
  function get_Ecoul(Z,v,npo,ne) result(Ec) 
    real(8), intent(in) :: Z,v,npo,ne 
    real(8), dimension(3) :: Ec 
    real(8) :: rnuc,u,ufunc,dufuncdlnu 
    
    !Vws = (Z - v*npo) 
    rnuc = (3.d0*v/(4.d0*PI))**(1.d0/3.d0) 
    u = (ne - npo)/(Z/v - npo) 
    ufunc = 1.d0 - 1.5d0*u**(1.d0/3.d0) + 0.5d0*u 
    dufuncdlnu = 0.5d0*(u - u**(1.d0/3.d0))
    
    Ec(1) = 3.d0/5.d0/137.d0/rnuc*(Z-v*npo)**2
    Ec(2) = Ec(1)/v*(dufuncdlnu*Z/(Z-npo*v)  - ufunc/3.d0 - 2.d0*ufunc*npo*v/(Z-npo*v)) 
    Ec(3) =-Ec(1)*dufuncdlnu*ne/(ne-npo + 1.d-30) ! There is something wrong here
    Ec(1) = Ec(1)*ufunc 

  end function get_Ecoul

end module eos_nse_mod 

