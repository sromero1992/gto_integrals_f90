subroutine get_energy(ener)
  use com_vars
  implicit none
  real(8), intent(inout)      :: ener
  real(8)                     :: ener2
  !real(8)                            
  !D,H0,F
  ener=0.0d0
  ener2=0.0d0
  do i = 1, nbas
     do j = 1, nbas
        ener = ener + D_MAT(i,j)*( H0_MAT(i,j) + F_MAT(i,j) )
     end do
     ener2 = ener2 + D_MAT(i,i)*( H0_MAT(i,i) + F_MAT(i,i) )
  end do
  write(*,*) 'ENERGY IN GET_ENERGY: ', ener, ener2
end subroutine
