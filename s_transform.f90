subroutine s_transform(mode,verbose,X)
  use com_vars, only       : nbas, eig, eigV, eig_tmp, eigV_tmp, X_MAT
  real(8), intent(inout)  :: X(nbas,nbas)
  integer                 :: i, j, k, mode
  logical                 :: verbose
  !This subroutine is the transformation to get molecular orbitals
  ! mode 1 is for passing transformation to X_MAT=V*D*V'
  ! mode 2 is for passing transformation to     X=V*D*V'

  !S transformation
  if (mode .EQ. 1) then
     allocate(X_MAT(nbas,nbas))
     X_MAT(:,:) = 0.0d0
     do i = 1, nbas
        do j = 1, nbas
           do k = 1, nbas
              X_MAT(i,j) = X_MAT(i,j) + eigV(i,k)*eig(k)**(-0.5d0)*eigV(j,k)
           end do
        end do
     end do
     if (verbose) then
        write(*,*)
        write(*,*) '#################### Writting X_MAT ############################'
        do i = 1, nbas
          write(*,*) X_MAT(i,:)
        end do
     end if
     X=X_MAT
  else
     X(:,:) = 0.0d0   
     do i = 1, nbas
        do j = 1, nbas
           do k = 1, nbas
              !To get tmp, run dens for mode 2 to get tmp arrays
              X(i,j) = X(i,j) + eigV_tmp(i,k)*eig_tmp(k)**(-0.5d0)*eigV_tmp(j,k)
           end do
        end do
     end do
     if (verbose) then 
        write(*,*)
        write(*,*) '##################### Writting X_TMP  ###########################'
        do i = 1, nbas
          write(*,*) X(i,:)
        end do
     end if 
  end if 
end  subroutine
