program test_fermi

  implicit none 
  real(8) :: eta, ife, ufe, lfe, delta, dfe, dfea, max_err, err, f12, fm12
  double precision :: ifermi12, zfermi32, zfermi12, zfermim12
  integer :: i  
  
  eta = 1.d0
  delta = 1.d-5
  max_err = 0.d0
   
  do i=-100, 100 
    eta = i/1.d3 + 6.d0 
    dfe = ifermi12(eta*(1.d0 + delta), dfea) - ifermi12(eta*(1.d0-delta), dfea) 
    dfe = dfe/(2.d0*delta*eta)
    ife = ifermi12(eta, dfea)
    f12 = zfermi12(ife)
    fm12 = zfermim12(ife) 
    err = abs((dfe-dfea)/dfea)
    max_err = max(max_err, err) 
    write(6,'(15es10.3)')eta, ife, err, dfea, dfe, (2.0/fm12 - dfea)/dfea
  enddo
  
  if (max_err>1.d-9) stop 1
  
  stop 0 
end program test_fermi
