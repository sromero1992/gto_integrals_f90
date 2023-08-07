subroutine  dens(mode, verbose, E_VEC, D_TMP) 
!mode = 1 stores density in D_MAT
!mode = 2 stores density in D_TMP
  use module_com
  use diag 
  implicit none
  real(8),intent(inout) :: D_TMP(nbas,nbas), E_VEC(nbas,nbas)
  integer               :: mode
  logical               :: verbose

  write(*,*) 'BUILDING DENSITY'  
  if (verbose) then  
     write(*,*)
     write(*,*)'EVECTORS FOR DENSITY (COLUMN WISE) :'
     do i = 1, nbas
        write(*,*)E_VEC(i,:)
     end do 
  end if 

  nelec= int( el(1) + el(2) )
  write(*,*) 'TOTAL NUMBER OF ELECTRONS : ',nelec

  if (mode .EQ. 1) then
     D_MAT(:,:)=0.0d0
     open(15,file='DENSITY.txt')
     do i = 1, nbas
        do j = 1, nbas
           do k = 1, nelec/2
              D_MAT(i,j) = D_MAT(i,j) + E_VEC(i,k)*E_VEC(j,k) 
           end do
        end do
        write(15,*) D_MAT(i,:)
     end do
     close(15)
     if (verbose) then  
        write(*,*)
        write(*,*)'DENSITY MATRIX (D_MAT) :'
        do i = 1, nbas
           write(*,*)D_MAT(i,:)
        end do 
     end if 
  else
     D_TMP(:,:)=0.0d0
     do i = 1, nbas
        do j = 1, nbas
           do k = 1, nelec/2
              D_TMP(i,j) = D_TMP(i,j) + E_VEC(i,k)*E_VEC(j,k) 
           end do
        end do
     end do
  end if 

  write(*,*) 'END OF DENSITY BUILDING'
end subroutine dens

