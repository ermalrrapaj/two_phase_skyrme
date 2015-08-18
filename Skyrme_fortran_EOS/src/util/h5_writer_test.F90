program test_h5_writer
  
  use hdf5
  use h5_writer
  
  implicit none 
  
  integer :: i,j,k 
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
        
        pp(i,j,k) = DBLE(i*j*k) 
           
      enddo
    enddo 
  enddo
  
  
  call h5open_f(error)
  h5filename = trim(adjustl("eostable_test.h5"))
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
  ierr = ierr + h5_output_double(file_id,"pp",pp,dims3,rank) 
  
  call h5fclose_f(file_id,error)
  call h5close_f(error)
   
end program test_h5_writer
