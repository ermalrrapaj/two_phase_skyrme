module eos_gibbs_mod
  use eos_skyrme_mod
  use eos_com_mod 
  use nrtype 
  use nr 

  implicit none    
  
  type :: gibbs_region
    ! Holds data about a mixed phase region
    real(8) :: T 
    integer :: npoints
    real(8), dimension(:), allocatable :: mun,mup,pp
    real(8), dimension(:), allocatable :: npl,nnl,npu,nnu
  contains 
    procedure :: boundary => Gibbs_boundary 
    procedure :: get_point => get_two_phase_point
  end type gibbs_region

  private :: Gibbs_boundary 
  private :: get_two_phase_point

contains
  
  function get_two_phase_point(self, eos_in) result(eos_out) 
    ! Find the EOS of a point in the mixed phase (or not in the mixed phase)
    class(gibbs_region), intent(in) :: self 
    type(eos_com), intent(in) :: eos_in 
    type(eos_com) :: eos_out
    
    type(eos_com) :: eos,eos_lo,eos_hi
    real(8), dimension(self%npoints) :: nplam 
    real(8) :: u,npl,nnl,npu,nnu,ulam,del,nbl,mun,xg(3),ff(3)
    integer :: i,ilo,imid,ihi 
    logical :: bad 
    
    if (.not. eos_in%densities_set) then
      write(6,"(a)")"Densities not set on call to get_two_phase_point"
    endif 
     
    eos_out = eos_in
    
    ! Make a list of proton densities for chosen densities  
    do i=1,self%npoints 
      ulam = (eos_out%nb - (self%npl(i) + self%nnl(i)))/((self%npu(i) + self%nnu(i)) - (self%npl(i) + self%nnl(i)))
      nplam(i) = ulam*self%npu(i) + (1.d0-ulam)*self%npl(i) 
    enddo   
    
    ! Search for the desired proton fraction in our list 
    ihi = self%npoints
    ilo = 1
    imid = FLOOR(0.5d0*(ihi + ilo))  
    do i=1,50 
      if (nplam(imid)>=eos_out%np) then
        ilo = imid
        imid = FLOOR(0.5d0*(ihi + ilo))  
      else 
        ihi = imid
        imid = FLOOR(0.5d0*(ihi + ilo))  
      endif       
      if (ihi-ilo<2) exit  
    enddo   
    
    ! Get a guess for the properties of the mixed phase from our table 
    del = (eos_out%np - nplam(ilo))/(nplam(ihi)-nplam(ilo))  
    npl = self%npl(ihi)*del + self%npl(ilo)*(1.d0-del) 
    nnl = self%nnl(ihi)*del + self%nnl(ilo)*(1.d0-del) 
    npu = self%npu(ihi)*del + self%npu(ilo)*(1.d0-del) 
    nnu = self%nnu(ihi)*del + self%nnu(ilo)*(1.d0-del) 
    u = (eos_out%nb - (npl + nnl))/(npu + nnu - npl - nnl) 
    
    ! Determine if we should do a mixed phase search 
    ! The current bounds on u are arbitrary 
    if (u>-0.1d0 .and. u<1.1d0) then
      ! These are used to scale the NR equations 
      eos_out%pp  = self%pp(ihi)*del  + self%pp(ilo)*(1.d0-del) 
      eos_out%mun = self%mun(ihi)*del + self%mun(ilo)*(1.d0-del) 
      eos_out%mup = self%mup(ihi)*del + self%mup(ilo)*(1.d0-del) 
      eos_out%u   = u 
      ! Set our initial NR guess 
      xg(1) = LOG(npu) 
      xg(2) = LOG(npl) 
      xg(3) = eos_out%mun 
      ! Perform NR iteration to solve phase equilibrium equations
      call newt(xg,bad,get_f)
      if (bad) then 
        print *,'bad point'
        ff = get_f(xg) 
        print *,xg,ff
      endif   
      npu = EXP(xg(1)) 
      npl = EXP(xg(2)) 
      mun = xg(3)
      
      ! Calculate properties of two phase medium 
      eos_lo = get_bulk_state_from_mixed(eos_input_from_mixed(npl,mun,eos_out%T))
      eos_hi = get_bulk_state_from_mixed(eos_input_from_mixed(npu,mun,eos_out%T))
      eos_out%u = (eos_out%nb - eos_lo%nb)/(eos_hi%nb - eos_lo%nb) 
      if (eos_out%u<0.d0 .or. eos_out%u>1.d0) then
        ! We are actually outside of the two phase region so just do single phase
        eos_out = get_bulk_state(eos_out)
        eos_out%npl = eos_out%np
        eos_out%nnl = eos_out%nn
        eos_out%npu = eos_out%np
        eos_out%nnu = eos_out%nn
        eos_out%u   = MIN(1.d0,MAX(u,0.d0))
      else   
        ! We are in the two phase region  
        eos_out%mixed_phase = .true.
        eos_out%pp  = eos_lo%pp  
        eos_out%mup = eos_lo%mup  
        eos_out%mun = eos_lo%mun  
        eos_out%mue = eos_lo%mue
        eos_out%ss  = (1.d0-eos_out%u)*eos_lo%nb*eos_lo%ss/eos_out%nb
        eos_out%ss  = eos_out%ss + eos_out%u*eos_hi%nb*eos_hi%ss/eos_out%nb 
        eos_out%ee  = (1.d0-eos_out%u)*eos_lo%nb*eos_lo%ee/eos_out%nb
        eos_out%ee  = eos_out%ee + eos_out%u*eos_hi%nb*eos_hi%ee/eos_out%nb 
        eos_out%npl = eos_lo%np
        eos_out%nnl = eos_lo%nn
        eos_out%npu = eos_hi%np
        eos_out%nnu = eos_hi%nn
      endif 
    else  
      eos_out = get_bulk_state(eos_out)
      eos_out%mixed_phase = .false. 
      eos_out%npl = eos_out%np
      eos_out%nnl = eos_out%nn
      eos_out%npu = eos_out%np
      eos_out%nnu = eos_out%nn
      eos_out%u   = MIN(1.d0,MAX(u,0.d0))
    endif 
    
  contains 
    
    function get_f(x_,output_) result(y_)
      ! System of equations to zero to find phase equilibrium for 
      ! fixed densities and temperature
      real(SP), intent(IN), dimension(:) :: x_
      logical, optional :: output_
      real(SP), dimension(size(x_)) :: y_
      type(eos_com) :: eosL_,eosH_
      real(8) :: u_,npl_,npu_,mun_
      npu_ = EXP(x_(1))
      npl_ = EXP(x_(2))
      mun_ = x_(3) 
      u_ = (eos_out%np - npl_)/(npu_ - npl_) 
      eosL_ = get_bulk_state_from_mixed(eos_input_from_mixed(npl_,mun_,eos_out%T))
      eosH_ = get_bulk_state_from_mixed(eos_input_from_mixed(npu_,mun_,eos_out%T)) 
      y_(1) = (eosL_%pp  - eosH_%pp )/eos_out%pp
      y_(2) = (eosL_%mup - eosH_%mup)/eos_out%mup
      y_(3) = (u_*eosH_%nn +  (1.d0-u)*eosL_%nn)/eos_out%nn - 1.d0 
    end function get_f 

  end function get_two_phase_point

  subroutine Gibbs_boundary(self,T)
    !  Find Gibbs phase boundaries for a fixed temperature
    class(gibbs_region), intent(OUT) :: self  
    real(8), intent(IN) :: T
    
    real(8), dimension(1000) :: munt,mupt,ppt,npl,npu,nnl,nnu  
    real(8) :: mun,xgt(2),xgo(2),fac
    integer :: i,j,idx,idxT,idxmax 
    logical :: worked 
    character(80) :: fname 
  
    ! Iterate over initial grid 
    idx = 100
    do j=1,801 
      mun = (DBLE(j-1)/800.d0*100.d0 - 50.1d0)*MAX(1.d0/197.3d0,0.d0)
      call try_to_find_solution(mun,xgo,1,xgt,worked)  
      if (worked) xgo = xgt    
    enddo 
     
    ! Work up to higher chemical potential
    xgo(1) = LOG(npl(idx)) 
    xgo(2) = LOG(npu(idx))
    fac = 3.d-2 
    do i=1,100 
      if (munt(idx)>0.d0) then
        mun = munt(idx)*(1.d0 + fac)
      else 
        mun = munt(idx)*(1.d0 - fac)
      endif 
      call try_to_find_solution(mun,xgo,1,xgt,worked)  
      if (worked) then
        xgo = xgt   
      else 
        fac = fac*0.5d0
      endif  
    enddo
    idxmax = idx 
     
    ! Try to take some small steps down
    idx = 101 
    fac = 5.d-2 
    do i=1,200 
      if (munt(idx)<0.d0) then
        mun = munt(idx)*(1.d0 + fac)
      else 
        mun = munt(idx)*(1.d0 - fac)
      endif 
      
      call try_to_find_solution(mun,xgo,-1,xgt,worked)  
      if (idx==1) exit 
      if (worked) then
        xgo = xgt   
      else 
        fac = fac*0.5d0
      endif  
    enddo
    
    ! Save the phase boundaries to type arrays 
    self%npoints = idxmax-idx+1
    if (ALLOCATED(self%pp)) DEALLOCATE(self%pp,self%mun,self%mup,self%nnl,self%nnu,self%npl,self%npu)  
    ALLOCATE(self%pp(self%npoints),self%mun(self%npoints),self%mup(self%npoints),self%nnl(self%npoints))
    ALLOCATE(self%nnu(self%npoints),self%npl(self%npoints),self%npu(self%npoints))  
    self%T = T
    self%pp(1:self%npoints) = ppt(idx:idxmax) 
    self%mun(1:self%npoints) = munt(idx:idxmax) 
    self%mup(1:self%npoints) = mupt(idx:idxmax) 
    self%nnl(1:self%npoints) = nnl(idx:idxmax) 
    self%nnu(1:self%npoints) = nnu(idx:idxmax) 
    self%npl(1:self%npoints) = npl(idx:idxmax) 
    self%npu(1:self%npoints) = npu(idx:idxmax) 
    
    return 
  
  contains
  
    subroutine try_to_find_solution(mun,xgin,increment,xgout,success)  
      real(8), intent(in) :: mun,xgin(2)
      integer, intent(in) :: increment
      real(8), intent(out) :: xgout(2)
      logical, intent(out) :: success  
      real(8) :: xg(2) 
      logical :: bad 
      type(eos_com) :: eos(2)
       
      success = .false. 
  
      ! First find the zero pressure point at high density 
      xg(2) = LOG(0.2d0) 
      call newt(xg(2:2),bad,get_zero_P)
      
      ! Now find the low density point with the same chemical potentials
      eos(1) = get_bulk_state_from_mixed(eos_input_from_mixed(EXP(xg(2)),mun,T))    
      eos = get_bulk_state_from_mu(eos_input_from_mu(eos(1)%mup,eos(1)%mun,T))    
      xg(1) = LOG(eos(1)%np) 
      
      ! Use these guesses to solve
      call newt(xg,bad,get_f)
      eos(1) = get_bulk_state_from_mixed(eos_input_from_mixed(EXP(xg(1)),mun,T))    
      eos(2) = get_bulk_state_from_mixed(eos_input_from_mixed(EXP(xg(2)),mun,T))    
      
      ! Ok, that didn't work, try the last succesful point if it exists
      if ((bad .or. ABS(eos(1)%np/eos(2)%np - 1.d0)<1.d-5) ) then
        xg = xgin
        call newt(xg,bad,get_f)
        eos(1) = get_bulk_state_from_mixed(eos_input_from_mixed(EXP(xg(1)),mun,T))    
        eos(2) = get_bulk_state_from_mixed(eos_input_from_mixed(EXP(xg(2)),mun,T))    
      endif
  
      if (ABS(eos(1)%nb/eos(2)%nb - 1.d0)>1.d-6 .and. ABS(eos(1)%xp/eos(2)%xp-1.d0)>1.d-6 & 
          .and. ABS(eos(1)%mup/eos(2)%mup - 1.d0)<1.d-6 & 
          .and. ABS(eos(1)%pp/eos(2)%pp - 1.d0)<1.d-6 &
          .and. eos(1)%nb < eos(2)%nb & 
          .and. eos(1)%full_set .and. eos(2)%full_set ) then
        xgout = xg
        success = .true.  
        write(6,'(20E16.3e3)',ADVANCE='no')mun/T,eos(1)%np,eos(1)%nn,eos(2)%np,eos(2)%nn,eos(1)%pp,eos(1)%mup/T
        write(6,'(i4)')idx
        
        idx = idx+increment
        
        if (idx<=SIZE(munt)) then 
          munt(idx) = mun 
          mupt(idx) = eos(1)%mup 
          ppt(idx)  = eos(1)%pp 
          npl(idx)  = eos(1)%np
          npu(idx)  = eos(2)%np
          nnl(idx)  = eos(1)%nn
          nnu(idx)  = eos(2)%nn
        endif 
      endif
    
    end subroutine try_to_find_solution 
  
    function get_zero_P(x_,output_) result(y_)
      real(SP), intent(IN), dimension(:) :: x_
      logical, optional :: output_
      real(SP), dimension(size(x_)) :: y_
      type(eos_com) :: eos_
      
      eos_ = get_bulk_state_from_mixed(eos_input_from_mixed(EXP(x_(1)),mun,T))
      y_(1) = eos_%pp/(eos_%nb*eos_%T) 
    
    end function get_zero_P
     
    function get_f(x_,output_) result(y_)
      real(SP), intent(IN), dimension(:) :: x_
      logical, optional :: output_
      real(SP), dimension(size(x_)) :: y_
      type(eos_com) :: eosL_,eosH_
       
      eosL_ = get_bulk_state_from_mixed(eos_input_from_mixed(EXP(x_(1)),mun,T))
      eosH_ = get_bulk_state_from_mixed(eos_input_from_mixed(EXP(x_(2)),mun,T))
      
      y_(1) = (eosL_%pp  - eosH_%pp )/(ABS(eosL_%pp)  + ABS(eosH_%pp))
      y_(2) = (eosL_%mup - eosH_%mup)/(ABS(eosL_%mup) + ABS(eosH_%mup)) 
    
    end function get_f 
  
  end subroutine gibbs_boundary
end module eos_gibbs_mod
