program test_nse
  
  use eos_skyrme_mod
  use eos_com_mod
  use hdf5
  use h5_writer
  implicit none 
  integer :: i,j,k 
  type(eos_com) :: eos 
  real(8) :: T,np,nn,etan

  integer, parameter :: nT  = 10
  integer, parameter :: nnp = 250
  integer, parameter :: nnn = 250 
  
  real(8), dimension(nnn)  :: nn_arr,etan_arr 
  real(8), dimension(nnp)  :: np_arr 
  real(8), dimension(nT)  :: T_arr 
  real(8), dimension(nnp,nnp,nT) :: mun,mup,pp,u
  real(8), dimension(nnp,nnp,nT,7) :: xx
  
  integer error, rank,ierr
  integer(HID_T) file_id,dset_id,dspace_id,aspace_id,attr_id
  integer(HID_T) atype_id
  integer(SIZE_T) attrlen

  integer(HSIZE_T) dims1(1), dims3(3), dims4(4), dims6(6)

  character(len=128) h5filename
  !T = 1.d0/197.3d0 
  !np = 0.9d-3
  !nn = 1.0d-3  
  !eos = get_bulk_state(eos_input_from_mixed(np,-T*0.1d0,T),.true.) 
  !return 
  do i=1,nnp
    do j=1,nnn 
      do k=nT,1,-1 
        
        np   = 10.d0**(DBLE(i-1)/DBLE(nnp-1)*(LOG10(0.4d0) - LOG10(1.d-10)) + LOG10(1.d-10)) 
        nn   = 10.d0**(DBLE(j-1)/DBLE(nnn-1)*(LOG10(1.d0) - LOG10(1.d-10)) + LOG10(1.d-10)) 
        etan = (DBLE(j-1)/DBLE(nnn-1)*100.d0 - 50.d0) 
        T  = 10.d0**(DBLE(k-1)/DBLE(nT-1)*(LOG10(20.d0) - LOG10(1.d0)) + LOG10(1.d0))/197.3d0  
        
        np_arr(i)  = np
        nn_arr(j)  = nn
        etan_arr(j) = etan  
        T_arr(k)   = T

        ! Right now the EoS is perfectly symmetric
        if (j>=i) then
          !eos = get_bulk_state(eos_input_from_densities(np,nn,T)) 
          eos = get_bulk_state(eos_input_from_mixed(np,etan*T,T)) 
          if (.not. eos%full_set) then 
            write(6,'(a,20E12.3)')'Bad point',np,nn,T*197.3d0
            write(6,'(30E12.3)')np,eos%np,eos%nn,T*197.3d0,eos%xx(1:7),SUM(eos%xx(1:7))
          endif 
          
          mun(i,j,k) = eos%mun 
          mup(i,j,k) = eos%mup
          pp(i,j,k)  = eos%pp
          u(i,j,k)   = eos%u
          xx(i,j,k,1:7) = eos%xx(1:7)
           
          mun(j,i,k) = eos%mup
          mup(j,i,k) = eos%mun
          pp(j,i,k)  = eos%pp
          u(j,i,k)   = eos%u
          xx(j,i,k,1:7) = eos%xx(1:7) 
        endif  
         
        !if (100==j.and.k==1) write(6,'(30E12.3)')np,eos%np,eos%nn,T*197.3d0,eos%mun,eos%mup,eos%pp,eos%xx(1:7),SUM(eos%xx(1:7))
        if (100==j.and.k==1) write(6,'(30E12.3)')np,eos%np,eos%nn,T*197.3d0,etan,eos%mun/T,eos%mup/T,eos%pp
           
      enddo
    enddo 
  enddo
  
  
  call h5open_f(error)
  h5filename = trim(adjustl("eostable_Skyrme.h5"))
  call h5fcreate_f(h5filename,H5F_ACC_TRUNC_F,file_id,error)
  
  rank=1
  dims1(1) = nnp
  ierr = ierr + h5_output_double(file_id,"np",np_arr,dims1,rank) 
  rank=1
  dims1(1) = nnn
  ierr = ierr + h5_output_double(file_id,"nn",nn_arr,dims1,rank) 
  rank=1
  dims1(1) = nnn
  ierr = ierr + h5_output_double(file_id,"etan",etan_arr,dims1,rank) 
  rank=1
  dims1(1) = nT
  ierr = ierr + h5_output_double(file_id,"T",T_arr,dims1,rank)
  
  rank=3
  dims3(1) = nnp
  dims3(2) = nnn
  dims3(3) = nT
  ierr = 0
  ierr = ierr + h5_output_double(file_id,"mup",mup,dims3,rank) 
  ierr = ierr + h5_output_double(file_id,"mun",mun,dims3,rank) 
  ierr = ierr + h5_output_double(file_id,"pp",pp,dims3,rank) 
  ierr = ierr + h5_output_double(file_id,"u",u,dims3,rank) 
  
  rank=4
  dims4(1) = nnp
  dims4(2) = nnn
  dims4(3) = nT
  dims4(4) = 7
  ierr = ierr + h5_output_double(file_id,"xx",xx,dims4,rank) 
   
  call h5fclose_f(file_id,error)
  call h5close_f(error)
   
end program test_nse 
