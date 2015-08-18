module tabular_EOS
  use eos_com_mod
  use HDF5
  use H5LT
   
  implicit none 
  
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
  
  ! Table data
  logical, private :: table_loaded = .false.  
  real(8), dimension(:,:,:), allocatable, private :: mun_arr,mup_arr,pp_arr
  real(8), dimension(:), allocatable :: np_arr,nn_arr,T_arr
  real(8), private :: dlnp,dlnn,dlT  
  real(8), private :: lnp0,lnn0,lT0
    
contains
  subroutine read_table(fname) 
    
    character(*), intent(in) :: fname 

    character(LEN=80) :: dset_name
    integer           :: error,rank
    integer(HID_T)    :: fid,did
    integer(HSIZE_T)  :: len,dims1(1),dims2(2),dims3(3),dims4(4),dims6(6)
    integer(SIZE_T)   :: type_size     
    integer           :: type_class
    
    call h5open_f(error)
    call h5fopen_f(TRIM(fname),H5F_ACC_RDONLY_F,fid,error)
    
    call h5ltget_dataset_info_f(fid,'np',dims1, type_class, &
                                type_size, error)
    allocate(np_arr(dims1(1))) 
    call h5ltread_dataset_double_f(fid,'np',np_arr,dims1,error)
    
    call h5ltget_dataset_info_f(fid,'nn',dims1, type_class, &
                                type_size, error)
    allocate(nn_arr(dims1(1))) 
    call h5ltread_dataset_double_f(fid,'nn',nn_arr,dims1,error)
     
    call h5ltget_dataset_info_f(fid,'T',dims1, type_class, &
                                type_size, error)
    allocate(T_arr(dims1(1))) 
    call h5ltread_dataset_double_f(fid,'T',T_arr,dims1,error)
     
    call h5ltget_dataset_info_f(fid,'pp',dims3, type_class, &
                                type_size, error)
    allocate(mup_arr(dims3(1),dims3(2),dims3(3)))  
    allocate(mun_arr(dims3(1),dims3(2),dims3(3)))  
    allocate(pp_arr(dims3(1),dims3(2),dims3(3)))  
    
    call h5ltread_dataset_double_f(fid,'mup',mup_arr,dims3,error)
    call h5ltread_dataset_double_f(fid,'mun',mun_arr,dims3,error)
    call h5ltread_dataset_double_f(fid,'pp',pp_arr,dims3,error)
    
    call h5fclose_f(fid,error)
    call h5close_f(error)
     
    dlnp = LOG(np_arr(2)) - LOG(np_arr(1))  
    dlnn = LOG(nn_arr(2)) - LOG(nn_arr(1))  
    dlT  = LOG(T_arr(2))  - LOG(T_arr(1))  
    lnp0 = LOG(np_arr(1)) 
    lnn0 = LOG(nn_arr(1)) 
    lT0  = LOG(T_arr(1)) 

    table_loaded = .true.
     
    return 
  end subroutine read_table 
   
  subroutine get_zero_bulk_mu(nt,xp) 
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
    type(eos_com), intent(in) :: eos_in 
    type(eos_com) :: eos_out,eos_uni
    logical :: bad,output
    real(8) :: lnp,lnn,lT,xnp,xnn,xT 
    integer :: inp,inn,iT 
     
    if (.not. table_loaded) then
      print *,'You havent loaded an EoS' 
      stop
    endif  
    
    eos_out = eos_in 
    
    lnp = LOG(eos_in%np) 
    lnn = LOG(eos_in%nn) 
    lT  = LOG(eos_in%T) 
 
    inp = FLOOR((lnp - lnp0)/dlnp)
    inn = FLOOR((lnn - lnn0)/dlnn)
    iT  = FLOOR((lT  - lT0)/dlT) 
    if (inp<1) inp = 1 
    if (inn<1) inn = 1 
    if (iT<1)  iT  = 1 
    if (inp>SIZE(np_arr)-1) inp = SIZE(np_arr)-1
    if (inn>SIZE(nn_arr)-1) inn = SIZE(nn_arr)-1
    if (iT> SIZE(T_arr)-1)  iT  = SIZE(T_arr)-1
    
    xnp = (lnp - LOG(np_arr(inp)))/dlnp 
    xnn = (lnn - LOG(nn_arr(inn)))/dlnn 
    xT  = (lT - LOG(T_arr(iT)))/dlT 
    
    eos_out%mup = threed_interp((/xnp,xnn,xT/),(/inp,inn,iT/),mup_arr) 
    eos_out%mun = threed_interp((/xnp,xnn,xT/),(/inp,inn,iT/),mun_arr) 
    eos_out%pp  = threed_interp((/xnp,xnn,xT/),(/inp,inn,iT/),pp_arr) 
    eos_out%full_set = .true.
    
  contains 
    function threed_interp(x,idx,arr) result(yy) 
      real(8), intent(in) :: x(:),arr(:,:,:) 
      integer, intent(in) :: idx(:) 
      real(8) :: yy 
      
      yy = arr(idx(1)  ,idx(2)  ,idx(3)  )*(1.d0-x(1))*(1.d0-x(2))*(1.d0-x(3)) & 
         + arr(idx(1)+1,idx(2)  ,idx(3)  )*(     x(1))*(1.d0-x(2))*(1.d0-x(3)) & 
         + arr(idx(1)  ,idx(2)+1,idx(3)  )*(1.d0-x(1))*(     x(2))*(1.d0-x(3)) & 
         + arr(idx(1)  ,idx(2)  ,idx(3)+1)*(1.d0-x(1))*(1.d0-x(2))*(     x(3)) & 
         + arr(idx(1)  ,idx(2)+1,idx(3)+1)*(1.d0-x(1))*(     x(2))*(     x(3)) & 
         + arr(idx(1)+1,idx(2)  ,idx(3)+1)*(     x(1))*(1.d0-x(2))*(     x(3)) & 
         + arr(idx(1)+1,idx(2)+1,idx(3)  )*(     x(1))*(     x(2))*(1.d0-x(3)) & 
         + arr(idx(1)+1,idx(2)+1,idx(3)+1)*(     x(1))*(     x(2))*(     x(3)) 
    
    end function threed_interp 
  end function get_bulk_state
  
  !function get_equilibrium_xp(nt) result(xp) 
  !  real(8) :: xp 
  !  real(8), intent(in) :: nt 
  !  real(8) :: ee,pp,muu,mul,muc,xu,xl 
  !  integer :: i 
  !  
  !  xu = 1.0d0 
  !  xl = 0.0d0 
  !  call get_state(nt,xu,ee,pp,muu,.true.)  
  !  call get_state(nt,xl,ee,pp,mul,.true.) 
  !  if (muu*mul>0.d0) then
  !    print *,'bad start' 
  !  endif 
  !  do i=1,50
  !    call get_state(nt,(xu+xl)*0.5d0,ee,pp,muc,.true.)  
  !    if (muc*muu>0.d0) then
  !      muu = muc 
  !      xu = 0.5d0*(xu+xl)  
  !    else 
  !      mul = muc 
  !      xl = 0.5d0*(xu+xl)  
  !    endif
  !  enddo
  !  if (abs(muc)>1.d-8) print *,'Search failed',muc 
  !  xp = xl     
  !end function
   
end module tabular_EOS 

