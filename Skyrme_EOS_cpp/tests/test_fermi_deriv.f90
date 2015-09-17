program test_fermi

  implicit none 
  real(8) :: eta, ife, ufe, lfe, delta, dfe, dfea, max_err, err
  double precision :: ifermi12, zfermi32 
  integer :: i  
  eta = 1.d0
  delta = 1.d-5
  max_err = 0.d0 
  do i=-100, 100 
    eta = i/1.d3 + 6.d0 
    dfe = zfermi32(eta*(1.d0 + delta), dfea) - zfermi32(eta*(1.d0-delta), dfea) 
    dfe = dfe/(2.d0*delta*eta)
    ife = zfermi32(eta, dfea)
    err = abs((dfe-dfea)/dfea)
    max_err = max(max_err, err) 
    print *,eta, ife, err, dfea, dfe
  enddo
  if (max_err>1.d-9) stop 1
  stop 0 
end program test_fermi
