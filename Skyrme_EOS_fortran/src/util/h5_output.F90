module h5_writer 
contains
  function h5_output_double(lfid,var_name,var_array,ldims,lrank) result(lierr)
    use hdf5
    implicit none
    
    integer :: lierr
    integer(HID_T), intent(in) :: lfid
    character(LEN=*), intent(in)          :: var_name
    real(8), dimension(*), intent(in) :: var_array
    integer(HSIZE_T), dimension(:),     intent(in) :: ldims
    integer, intent(in) :: lrank 

    integer(HID_T) :: ldspace_id,ldset_id,mem_type
    integer :: lerror 
    mem_type = H5T_NATIVE_DOUBLE
    call h5screate_simple_f(lrank, ldims, ldspace_id, lerror)
    lierr = lerror
    call h5dcreate_f(lfid, TRIM(var_name), mem_type, &
         ldspace_id, ldset_id, lerror)
    lierr = lerror + lierr
    call h5dwrite_f(ldset_id, mem_type, var_array, &
         ldims, lerror)
    lierr = lerror + lierr
    call h5dclose_f(ldset_id, lerror)
    lierr = lerror + lierr
    call h5sclose_f(ldspace_id, lerror)
    lierr = lerror + lierr
    if (lierr .ne. 0) then
      write(6,'(a,a,a)') 'Error writing variable ',TRIM(var_name),' to file.' 
    else
      write(6,'(a,a,a)') TRIM(var_name),' written to file.' 
    endif 
        
    return 
  end function h5_output_double

  function h5_output_int(lfid,var_name,var_array,ldims,lrank) result(lierr)
    use hdf5
    implicit none
    
    integer :: lierr
    integer(HID_T), intent(in) :: lfid
    character(LEN=*), intent(in)          :: var_name
    integer, dimension(*), intent(in) :: var_array
    integer(HSIZE_T), dimension(:),     intent(in) :: ldims
    integer, intent(in) :: lrank 

    integer(HID_T) :: ldspace_id,ldset_id,mem_type
    integer :: lerror 

    mem_type = H5T_NATIVE_INTEGER
    call h5screate_simple_f(lrank, ldims, ldspace_id, lerror)
    lierr = lerror
    call h5dcreate_f(lfid, TRIM(var_name), mem_type, &
         ldspace_id, ldset_id, lerror)
    lierr = lerror + lierr
    call h5dwrite_f(ldset_id, mem_type, var_array, &
         ldims, lerror)
    lierr = lerror + lierr
    call h5dclose_f(ldset_id, lerror)
    lierr = lerror + lierr
    call h5sclose_f(ldspace_id, lerror)
    lierr = lerror + lierr
    if (lierr .ne. 0) then
      write(6,'(a,a,a)') 'Error writing variable ',TRIM(var_name),' to file.' 
    else
      write(6,'(a,a,a)') TRIM(var_name),' written to file.' 
    endif 
        
    return 
  end function h5_output_int
end module h5_writer

!  ierr = ierr + h5_output_int(file_id,"num_isotopes",(/num_isotopes/),dims1,rank) 
!  
!  rank = 1
!  dims1(1) = nrho_nse
!  ierr = ierr + h5_output_double(file_id,"logrho",logrho_nse,dims1,rank) 
