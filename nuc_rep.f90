subroutine  nuc_rep(ener_nuc)
  use class_basis_function
  use mat_build
  use module_com
  use module_g
  implicit none
  real(8), intent(inout) :: ener_nuc 
  real(8)                :: d

  ener_nuc = 0.0d0
  if (all_atms .EQ. 1 ) then
     ener_nuc = 0.0d0
  else
     do i = 1, all_atms
        do j = i, all_atms
           if ( i  .NE. j) then
              d    = norm2( g(i)%origin(:) - g(j)%origin(:))
              ener_nuc = ener_nuc + 1.0d0*( g(i)%qn * g(j)%qn ) / d 
           end if
        end do
     end do
  end if


end subroutine
