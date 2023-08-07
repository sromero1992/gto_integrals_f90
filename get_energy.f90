subroutine get_energy(ener)
  use module_com
  implicit none
  real(8), intent(inout)      :: ener
  real(8)                     :: ener2,trace, H_tmp(nbas, nbas), DH_MAT(nbas,nbas)
  !real(8)                            
  !D,H0,F
  H_tmp = H0_MAT + F_MAT 
  call dgemm( 'N', 'N', nbas, nbas, nbas, 1.0d0, & 
              D_MAT, nbas, H_tmp, nbas, 0.0d0, DH_MAT, nbas)
  trace = 0.0d0
  do i = 1, nbas 
     trace = trace + DH_MAT(i,i)
  end do
  ener = trace
end subroutine
